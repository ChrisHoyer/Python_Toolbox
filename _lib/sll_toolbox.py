#############################################################################
#   Different Scripts to Model State-Locked-Loops
#
#   - Net_2Mutually: Simple Network in Laplace with two bidirectional coupled nodes to Calc Global Freq
#   - Perform_Timeseries_Simulation: Internal Simulation Engine
#
#   Prototype Models
#   - Calc_2Cplg_linear: Calculates important Values for Steady State Analysis
#   - Calc_linear: Calculates important Values for Steady State Analysis
#   - Calc_2Cplg_nonlinear: Calculates important Values for Steady State Analysis (nonlinear)
#   - Calc_nonlinear: Calculates important Values for Steady State Analysis (nonlinear)
#
#  Dynamical Models
#   - Run_Timeseries: Calulate Euler Simulation of Coupled Oscillators
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Date: 06-11-2021
#############################################################################
import controltheory_toolbox as ctrl
import basic_toolbox as basic

import sympy
import time
import numpy as np
import scipy as sp
from scipy import signal
import matplotlib.pyplot as plt

#############################################################################
###                 PLL Network Solver
#############################################################################

# Solving two mutually coupled PLL Nodes
def Net_2Mutually(phase_delay, omega0_div, G_CPLG, variable, 
                  mode = 0, nonlinear = False,
                  feedback_phase = 0,
                  coupling_function= "triangle",
                  Nyquist_Freq = np.linspace(0,1e6,int(1e5)),
                  K_LF = None,
                  Nyquist_Neglect = 10,
                  calc_networkstab = True,
                  nyquist_tolerance = 1e-3, Nyquist_Real = 0.999,
                  stab_tolerance = 1e-14,
                  debug_return = True, freqinHz = True, 
                  Export_OLTF = False):
#############################################################################    
    """
   Simple Network in Laplace with two bidirectional coupled nodes to calculate
   the global network frequency and stability vs the phase delay 
   Note, that this is a differential equation. Where the init value has to be set.
    
    parameters              description
    =====================  =============================================:
    phase_delay             current phase delay
    omega0_div              divided quiescent PLL frequency
    G_CPLG                  Open Loop Transfer Function in steady state (with laplace variable)
                            (nonlinear) if nonlinear selected: open loop transfer function in steady state (K_LF=1)
    variable                laplace variable
    
    nonlinear               (optionally) nonlinear coupling function?
    feedback_phase          (optionally) inverter in feedback path
    mode                    (optionally) in-phase or anti-phase locking mode
    coupling_function       (optionally) PD coupling function (XOR=sawtooth, Mixer=cos)
    Nyquist_Freq            (optionally) Span where Nyquist Stability is anaylized
    Nyquist_Neglect         (optionally) Neglect first Values for Analysis
    calc_networkstab        (optionally) Calculate Network Stability via Nyquist
    K_LF                    (optionally) LF Transfer function with laplace variable (no steady state)
    nyquist_tolerance       (optionally) tolerance of nyquist solving
    Nyquist_Real            (optionally) nyquist realpart checker
    stab_tolerance          (optionally) tolerance of stability solver
    debug_return            (optionally) Return optional values (see below)
    freqinHz                (optionally) Return all Frequencies in Hz

    return type: dict
    return parameters       description
    =====================  =============================================:       
    Network_Freq            Calculated Network Frequency in Radiants
    Network_Freq_Stable     Calculated stable Network Frequency in Radiants
    Network_Freq_UnStable   Calculated unstable Network Frequency in Radiants
    Time_Delay              Time Delay accoring to phase and Network_Freq
    Zeros                   Zeros for this Network Frequency Solution
    Poles                   Poles for this Network Frequency Solution
    Nyquist_Gain            Network Gain at Imaginary Part ~= 0
    
    H_PD                    (optionally) Used coupling function
    H_prime_PD              (optionally) derivertive of used coupling function
    Zeros_Success           (optionally) Bool for Solving Success of Zeros
    Poles_Success           (optionally) Bool for Solving Success of Poles    
    Nyquist_Solution        (optionally) Solution for Nyquist Stability
    Nyquist_Freq            (optionally) Corresponding Frequency Axis
    
    """

############################################################################# 


    # rescale return values to Hz
    if freqinHz:
        scaletoHz = 1/(2*np.pi)
    else:
        scaletoHz = 1
        
    # Use K_LF also for Network Stability
    if K_LF == None:
        K_LF = 1
    
    # return dictionary
    return_var = {}
    return_var["Phase_Delay"] = -phase_delay/np.pi
    return_var["Net_Calc"] = calc_networkstab

    # Calculate Phase Error tilde substitude
    phase_error = mode - phase_delay + feedback_phase       

#############################################################################   
#           Calculate Phase Error Transfer Function
############################################################################# 
    
    # coupling function, XOR
    if coupling_function == "triangle":
        
        # Calculate Triangle function        
        H_PD = 0.5*sp.signal.sawtooth(phase_error, width=0.5)
        
        # derivertive of the slope
        H_prime_PD = 1/np.pi*sp.signal.square(phase_error, duty=0.5)

    # coupling function, MIXER
    elif coupling_function == "-cos":
        
        # Calculate Mixer function        
        H_PD = -0.5 * np.cos(phase_error)
        
        # derivertive of the slope
        H_prime_PD = 0.5*np.sin(phase_error)
        
    # coupling function, PFD
    elif coupling_function == "sawtooth":
        
        # Calculate Sawthooth function        
        H_PD = 0.5 * sp.signal.sawtooth(phase_error)
        
        # derivertive of the slope
        H_prime_PD = 1/np.pi*sp.signal.square(phase_error, duty=1)  
        
    # return H_PD and H_prime_PD
    if debug_return:
        return_var["H_PD"] = H_PD
        return_var["H_prime_PD"] = H_prime_PD
        
#############################################################################   
#           Calculate Network Frequency
############################################################################# 

    if nonlinear:
        # initial value for solving for stability
        G_LS = G_CPLG(0)
        
        # calculate global frequency
        OmegaNet = G_CPLG(H_PD)
           
        # calculate time delay
        time_delay = phase_delay/OmegaNet   

# =========================================================================    
       
    else:
    
        # Locked State Gain (Steady State)
        G_LS = variable * G_CPLG
        G_LS = float(sympy.limit(G_LS, variable, 0))
        
        # calculate global frequency
        OmegaNet = G_LS * H_PD + omega0_div
        
        # calculate time delay
        time_delay = phase_delay/OmegaNet

# =========================================================================    
         
    # Return Values fpr global network frequency and time delay
    return_var["Network_Freq"] = OmegaNet*scaletoHz
    return_var["Time_Delay"] = time_delay
    
#############################################################################   
#           Calculate PLL Stability
#############################################################################  
 
    if nonlinear:
        # Transfer Function of Open Loop (derivative) with H_prime and H
        G_CPLG_TF = G_CPLG(H_PD, H_prime_PD=H_prime_PD, derivative=True) * 1/variable 
        
        # Transfer Function of Open Loop including LF
        G_CPLG_TF_LF = G_CPLG_TF * K_LF
        
 
# =========================================================================    
       
    else:
        
        # Transfer Function of Open Loop with H_prime
        G_CPLG_TF = G_CPLG * H_prime_PD
        
        # Transfer Function of Open Loop with LF with H_prime
        G_CPLG_TF_LF = G_CPLG_TF * K_LF
        
    
# =========================================================================    

    # Closed Loop Transfer Function of a Single PLL
    G_CPLG_CL = G_CPLG_TF / (1 + G_CPLG_TF)
    G_CPLG_CL = sympy.simplify(G_CPLG_CL)
    
    # Generate Nominator and Denmominator for Analyzes of Poles/Zeros
    G_CPLG_CL_nom, G_CPLG_CL_denom = sympy.fraction(G_CPLG_CL)   

    # export function for Denom and Nom
    G_CPLG_CL_nom = sympy.lambdify([variable], G_CPLG_CL_nom)
    G_CPLG_CL_denom = sympy.lambdify([variable], G_CPLG_CL_denom)
        
    # solving for stability
    solved_nom = sp.optimize.root(G_CPLG_CL_nom, G_LS, tol=stab_tolerance)
    solved_denom = sp.optimize.root(G_CPLG_CL_denom, G_LS, tol=stab_tolerance)

    # save solution
    zeros = solved_nom.x
    poles = solved_denom.x
 
    # check solution of zeros 
    nom_error = G_CPLG_CL_nom(solved_nom.x)
    if not(solved_nom.success) or not(nom_error < stab_tolerance):
        debug_info = "Phase: " + str(phase_delay) 
        debug_info = debug_info + "\t Nom Error: " + str(nom_error)   
        debug_info = debug_info + "\t Sol: " + str(solved_nom.x)             
        #print(debug_info)
          
   
    # check solution of poles 
    denom_error = G_CPLG_CL_denom(solved_denom.x)
    if not(solved_denom.success) or not(denom_error < stab_tolerance):
        debug_info = "Phase: " + str(phase_delay) 
        debug_info = debug_info + "\t Denom Error: " + str(denom_error) 
        debug_info = debug_info + "\t Sol: " + str(solved_denom.x) 
        #print(debug_info)
 
    # return zeros and poles
    return_var["Zeros"] = zeros    
    return_var["Poles"] = poles   

     # return Success and Stability
    if debug_return:
        
        # Check Pole for Stability
        return_var["Zeros_Success"] = solved_nom.success       
        return_var["Poles_Success"] = solved_denom.success  
         
        # Check Pole for Stability
        return_var["Stability"] = True if (poles < 0.0) else False
        
    # Return only Stable Frequencies
    return_var["Network_Freq_Stable"] = OmegaNet*scaletoHz if (poles < 0.0) else float("NaN")        
    return_var["Network_Freq_UnStable"] = OmegaNet*scaletoHz if (poles >= 0.0) else float("NaN") 

#############################################################################   
#           Calculate Loop Bandwidth
############################################################################# 
       
    if Export_OLTF:
                
        return_var["OpenLoopTF"] = G_CPLG_TF_LF
    
#############################################################################   
#           Calculate Network Stability
#############################################################################  
    
    # Export no Values:
    return_var["Nyquist_Point"] =  float("NaN") 
    return_var["Nyquist_Point_Freq"] =  float("NaN") 
    
    # return Nyquist Calculation
    return_var["Nyquist_Solution"] = float("NaN")
    return_var["Nyquist_Solution_Freq"] = float("NaN")

            
    # Calculate Network Stability
    if calc_networkstab and not(nonlinear):
        
        # Laplace Delay block for calculacted Delay
        exp_delay = sympy.exp(-time_delay*variable) 
            
        # Closed Loop Transfer Function of a Single PLL with LF
        G_CPLG_CL_LF = G_CPLG_TF_LF / (1 + G_CPLG_TF_LF)
        G_CPLG_CL_LF = sympy.simplify(G_CPLG_CL_LF)
        
        # Open loop network transfer function of coupled PLLs
        G_NET = G_CPLG_CL_LF * exp_delay  * G_CPLG_CL_LF * exp_delay
        
        # Calculate Nyquist Frequency Span
        Nyquist_Omega = 2j * np.pi * Nyquist_Freq
         
        try:
            # Nyquist Plot with Complex return
            Nyquist_Calc = ctrl.Extract_Sympy_1Var(G_NET, Nyquist_Omega,
                                                   variable=variable)
            # return Nyquist Calculation
            return_var["Nyquist_Solution"] = Nyquist_Calc
            return_var["Nyquist_Solution_Freq"] = Nyquist_Omega*scaletoHz
                                  
            # Find Indices where imaginary part sign is changed (except first)
            Nyquist_Sign = ((np.roll(np.sign(np.imag(Nyquist_Calc)), 1) - np.sign(np.imag(Nyquist_Calc))) != 0)
            Nyquist_Sign = [i for i, x in enumerate(Nyquist_Sign) if x and i > 0]
            
            # Only where real part is positive
            Nyquist_Sign = [i for i in Nyquist_Sign if np.real(Nyquist_Calc[i]) > Nyquist_Real]

            if len(Nyquist_Sign) > 1:
                Nyquist_Sign = [i for i in Nyquist_Sign if np.real(Nyquist_Calc[i]) > 0.5]
            
            if len(Nyquist_Sign) > 1:
                print("Error in Nyquist Point Calculation - Multiple Points close to 1+0j")
                    
            # Save Nyquist Point and Frequency from Real Part to 1
            if  len(Nyquist_Calc[Nyquist_Sign]) == 1:
                return_var["Nyquist_Point"] = Nyquist_Calc[Nyquist_Sign][0]       
                return_var["Nyquist_Point_Freq"] = np.imag(Nyquist_Omega[Nyquist_Sign][0] *scaletoHz)
            else:
                return_var["Nyquist_Point"] = float("NaN")       
                return_var["Nyquist_Point_Freq"]  = float("NaN")              
            
            
        except:
            print("Error in Nyquist Point Calculation")
            pass
            
            
    return return_var

# Iterative Task during Simulation
def Perform_Timeseries_Simulation(t_index, dt, results, Osc_Param, Network_Param,
                                  CheckEnableCplg=False):
    
    #############################################################################    
    """
    Subfunction, which is called for each transient step during simulation
    
    parameters              description
    =====================  =============================================:
    t_index                 current time index
    results                 result vector storing transient results
    Osc_Param               dictionary containing oscillator parameters
    Network_Param           dictionary containing network parameters
    CheckEnableCplg         bool to check if coupling can be enabled

    """
    #############################################################################
    
    def PD_Fct(phase1, phase2, osc_index):
        
        # phase error
        phase_error = phase1 - phase2
        
        # PD Gain
        pd_gain = Osc_Param[osc_index]["pd_amplitude"]
        
        if Osc_Param[osc_index]["pd_fct"].upper() == "XOR":
            return pd_gain * sp.signal.sawtooth(phase_error, width=0.5)
            
        if Osc_Param[osc_index]["pd_fct"].upper() == "SIN":
            return pd_gain * np.sin(phase_error)
           
        else:
            return pd_gain * np.sin(phase_error)

    #############################################################################
    
    def LF_Fct(V_out, V_in, osc_index):
      
        # first order RC
        if Osc_Param[osc_index]["lf_fct"].upper() == "1_RC":
            
            # PD Gain and time constant
            lf_gain = Osc_Param[osc_index]["lf_gain"]
            lf_tau = Osc_Param[osc_index]["lf_tau"]
            
            # calculate next step
            dVC1 = 1/(lf_tau[0]) * (V_in[t_index-1] - V_out[t_index-1])
            Vc1 = V_out[t_index-1] + dt * dVC1
            
            V_out[t_index] = lf_gain * Vc1
                      
            return V_out
        
        # second order RC
        if Osc_Param[osc_index]["lf_fct"].upper() == "2_RC":
            
            # additional parameters needed for second order
            if t_index == 1:
                results[osc_index]["lf_2nd"] = np.zeros(len(V_out))
            
            # PD Gain and time constant
            lf_gain = Osc_Param[osc_index]["lf_gain"]
            lf_tau = Osc_Param[osc_index]["lf_tau"]
            
            # calculate first filter
            dVC1 = 1/(lf_tau[0]) * (V_in[t_index-1] - results[osc_index]["lf_2nd"][t_index-1])
            results[osc_index]["lf_2nd"][t_index] = results[osc_index]["lf_2nd"][t_index-1] + dt * dVC1

            # calculate second filter
            dVC2 = 1/(lf_tau[1]) * (results[osc_index]["lf_2nd"][t_index] - V_out[t_index-1])
            Vc2 = V_out[t_index-1] + dt * dVC2
            
            V_out[t_index] = lf_gain * Vc2
                      
            return V_out        
        
        
        # no filter Vout = Vin
        else:
            return V_in

    #############################################################################
    
    def PD_Adder(osc_index):
        
        # calculate adder voltage
        v_add = 0
        
        # iterate each coupled node, do the PD and adding
        for cplg_index in range( np.size(Network_Param["cplg_gain"][osc_index, ::]) ):
            
            # Get coupling gain and delay
            cplg_gain = Network_Param["cplg_gain"][osc_index, cplg_index]
            cplg_delay = Network_Param["cplg_delay"][osc_index, cplg_index]
            
            # get latest values of coupled oscillator
            cplg_theta = results[cplg_index]["theta"]
            
            # calculate coupling phase value
            phase_own = results[osc_index]["theta"][ t_index - 1 ]
            phase_cplg = cplg_theta[ t_index - 1 - int(cplg_delay/dt) ]
            
            
            # check if coupling can be established, otherwise the input is zero
            if CheckEnableCplg and t_index <= int(cplg_delay/dt):
                    phase_cplg = 0
            
            # calculate coupling input
            v_add += cplg_gain * PD_Fct(phase_cplg, phase_own, osc_index)
            
        return v_add
    
    #############################################################################
    
    # Calculus for each oscillator
    for osc_index in range(len(results)):
          
        # adder voltage
        results[osc_index]["vadd"][t_index] = PD_Adder(osc_index)

        # calculate tuning voltage using loop filter
        results[osc_index]["vtune"] = LF_Fct( results[osc_index]["vtune"], results[osc_index]["vadd"], osc_index)
                
        # integration by the oscillator's 1/s
        osc_value = Osc_Param[osc_index]["omega0"] + Osc_Param[osc_index]["K_vco"] * results[osc_index]["vtune"][t_index]

        # Add oscillator noise to the phase
        if False:
             osc_value = osc_value + np.random.normal(0, 1e6)   
             
        # calulate instantaneous phase
        results[osc_index]["theta"][t_index] = results[osc_index]["theta"][t_index -1] + dt * osc_value
        
        # calculate frequency
        results[osc_index]["angular_frequency"][t_index] = (results[osc_index]["theta"][t_index] - results[osc_index]["theta"][t_index - 1])/dt
        
    
    # return results    
    return results
                            
#############################################################################
###              3th/3.5th Generation, Homogen and Linear
#############################################################################  
def Calc_2Cplg_linear(time_start=0, time_end=10e-9, time_points=1e3,
                      VCO_freq0 = 2*np.pi*24.25e9,
                      Div_N=16*32, K_VCO=2*np.pi*757.47e6,
                      K_PD = 1.6, K_LF = 1, fb_phase = 0,  
                      s = sympy.Symbol('s', complex=True),
                      Nqy_freq = np.linspace(1e3, 100e6, int(1e3)),
                      Nqy_skip = 10, Export_OLTF = False):
    #############################################################################    
    """
    Simple Network of two coupled PLL System in Laplace Domain with default values from
    3. Gen Prototypes (Version 1.1) with 2nd order RC-Loopfilter
    
    parameters              description
    =====================  =============================================:
    time_start              starting phase for iteration
    time_end                ending phase for iteration
    time_points             number of iterations
    
    VCO_freq0               VCO Closed Loop Freq (Hz) (default: 24.25e9)
    Div_N                   Divider Value (default: 16*32)
    K_VCO                   VCO Sensitivity (Hz/V) (default: 757e6)
    K_PD                    PD Sensitivity (V/pi) (default: 1.6)
    K_LF                    Loop Filter Function (default: 1, no laplace)
    s                       Laplace Variable (default: "s")
    Nqy_freq                Frequencies for which Laplace is calculated
    Nqy_skip                Time Delays for which the Nyquist Calulation is skipped
    
    fb_phase               (optionally) Feedback Phase
    """
    #############################################################################
    
    # Some Infos
    print("\033[34;1m" + "Simulation of two coupled Oscillators (without nonlinearity)")
    print("\033[34;1m" + "-> Sweep Time-delay between {}s to {}s with {} points".format(basic.EngNot(time_start), 
                                                                                        basic.EngNot(time_end), 
                                                                                        basic.EngNot(time_points)))
    print("\033[34;1m" + "-> Free running closed Loop Freq. {}Hz and Division by {}".format(basic.EngNot(VCO_freq0/(2*np.pi)),
                                                                                            Div_N))
    print("\033[34;1m" + "-> Sensitivity of VCO {}Hz/V and Phase-Detector {}V/pi".format(basic.EngNot(K_VCO/(2*np.pi)),
                                                                                         basic.EngNot(K_PD)))
    print("\033[34;1m" + "-> Perform Nyquist Calulations between {}Hz to {}Hz with {} points\033[0m".format(basic.EngNot(Nqy_freq[0]),
                                                                                                            basic.EngNot(Nqy_freq[-1]), 
                                                                                                            basic.EngNot(len(Nqy_freq))))
    
    start_time = time.time()
    
    #############################################################################
    #           Variables of PLL System and Network Analysis (Default Values)
    #############################################################################
    
    # Quiescent PLL Freq
    Omega0_Div = VCO_freq0/Div_N

    # Time Delay
    Time_Delay = np.linspace(time_start,time_end,int(time_points))
    
    # PD Phase-Error Transfer Function
    coupling_function = "triangle"
    
      
    #############################################################################
    #           Calculate PLL Solution and Network Stability
    #############################################################################
    
    # transfer function of PD Input to CrossCoupling output
    G_CPLG = K_PD * K_VCO/s * 1/Div_N
    
    # Mean Phase for Time Delay
    PhaseError = Omega0_Div * Time_Delay
             
    Solution_InPhase = []
    Solution_AntiPhase = []
    
    # Progressbar
    l = len(PhaseError)
    basic.printProgressBar(0, l, length = 50)
    
    for index, PhaseDelay in enumerate(PhaseError): 
       
        # calculate every 100th network stability
        calc_networkstab = True if (index%Nqy_skip == 0) else False
        
        # calculate in-phase
        Solution_InPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                               K_LF = K_LF,
                                               mode=0, coupling_function = coupling_function,
                                               Nyquist_Freq = Nqy_freq,
                                               calc_networkstab = calc_networkstab,
                                               feedback_phase = fb_phase,
                                               Export_OLTF = Export_OLTF))
        
        # calculate anti-phase    
        Solution_AntiPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                                 K_LF = K_LF,
                                                 mode=np.pi, coupling_function = coupling_function,
                                                 Nyquist_Freq = Nqy_freq,
                                                 calc_networkstab = calc_networkstab,
                                                 feedback_phase = fb_phase,
                                                 Export_OLTF = Export_OLTF))
        
       # Update Progress Bar
        basic.printProgressBar(index + 1, l, length = 50) 
         
    # Convert List with Dicts to Dicts with Lists
    Solution_InPhase = {k: [dic[k] for dic in Solution_InPhase] for k in Solution_InPhase[0]}
    Solution_AntiPhase = {k: [dic[k] for dic in Solution_AntiPhase] for k in Solution_AntiPhase[0]}
    
    # Some Infos
    elapsed_time = time.time() - start_time
    print("\033[0m" + "Calculations \033[32;1m" + "done" + "\033[0m! (Time needed: " + str(np.round(elapsed_time,3)) + "sec.)")
    
    # Generate VCO Frequency by Scaling with Division Factor
    Solution_InPhase["VCO_Freq_Stable"] = np.asarray(Solution_InPhase["Network_Freq_Stable"])*Div_N
    Solution_InPhase["VCO_Freq_UnStable"] = np.asarray(Solution_InPhase["Network_Freq_UnStable"])*Div_N   
    Solution_AntiPhase["VCO_Freq_Stable"] = np.asarray(Solution_AntiPhase["Network_Freq_Stable"])*Div_N
    Solution_AntiPhase["VCO_Freq_UnStable"] = np.asarray(Solution_AntiPhase["Network_Freq_UnStable"])*Div_N   
    
    # Some Infos
    print("\033[0m" + "Rescalings \033[32;1m" + "done" + "\033[0m!")
    
    #############################################################################    
    
    return {"m_0": Solution_InPhase,  "m_pi": Solution_AntiPhase}

def Calc_linear(time_start=0, time_end=10e-9, time_points=1e3,
                VCO_freq0 = 2*np.pi*24.25e9,
                Div_N=16*32, K_VCO=2*np.pi*757.47e6,
                K_PD = 1.6, K_LF = 1, fb_phase = 0,  
                s = sympy.Symbol('s', complex=True),
                Nqy_freq = np.linspace(1e3, 100e6, int(1e3)),
                Nqy_skip = 10, phi_mode = [0, np.pi], coupling_function = "triangle",
                Export_OLTF = False):
    #############################################################################    
    """
    Simple Network of two coupled PLL System in Laplace Domain with default values from
    3. Gen Prototypes (Version 1.1) with 2nd order RC-Loopfilter
    
    parameters              description
    =====================  =============================================:
    time_start              starting phase for iteration
    time_end                ending phase for iteration
    time_points             number of iterations
    
    VCO_freq0               VCO Closed Loop Freq (Hz) (default: 24.25e9)
    Div_N                   Divider Value (default: 16*32)
    K_VCO                   VCO Sensitivity (Hz/V) (default: 757e6)
    K_PD                    PD Sensitivity (V/pi) (default: 1.6)
    K_LF                    Loop Filter Function (default: 1, no laplace)
    s                       Laplace Variable (default: "s")
    Nqy_freq                Frequencies for which Laplace is calculated
    Nqy_skip                Time Delays for which the Nyquist Calulation is skipped
    
    fb_phase               (optionally) Feedback Phase
    """
    #############################################################################
    
    # Some Infos
    print("\033[34;1m" + "Simulation of two coupled Oscillators (without nonlinearity)")
    print("\033[34;1m" + "-> Sweep Time-delay between {}s to {}s with {} points".format(basic.EngNot(time_start), 
                                                                                        basic.EngNot(time_end), 
                                                                                        basic.EngNot(time_points)))
    print("\033[34;1m" + "-> Free running closed Loop Freq. {}Hz and Division by {}".format(basic.EngNot(VCO_freq0/(2*np.pi)),
                                                                                            Div_N))
    print("\033[34;1m" + "-> Sensitivity of VCO {}Hz/V and Phase-Detector {}V/pi".format(basic.EngNot(K_VCO/(2*np.pi)),
                                                                                         basic.EngNot(K_PD)))
    print("\033[34;1m" + "-> Perform Nyquist Calulations between {}Hz to {}Hz with {} points\033[0m".format(basic.EngNot(Nqy_freq[0]),
                                                                                                            basic.EngNot(Nqy_freq[-1]), 
                                                                                                            basic.EngNot(len(Nqy_freq))))
    
    start_time = time.time()
    
    #############################################################################
    #           Variables of PLL System and Network Analysis (Default Values)
    #############################################################################
    
    # Quiescent PLL Freq
    Omega0_Div = VCO_freq0/Div_N

    # Time Delay
    Time_Delay = np.linspace(time_start,time_end,int(time_points))
    
    # PD Phase-Error Transfer Function
    
      
    #############################################################################
    #           Calculate PLL Solution and Network Stability
    #############################################################################
    
    # transfer function of PD Input to CrossCoupling output
    G_CPLG = K_PD * K_VCO/s * 1/Div_N
    
    # Print Loop gain
    GLoop = float(s * G_CPLG)
    print("\033[34;1m" + "-> loop-gain: {}\033[0m".format(basic.EngNot(GLoop) ))                                          
                                                                                                        
    # Mean Phase for Time Delay
    PhaseError = Omega0_Div * Time_Delay
             
    Solution = {}
    
    # Progressbar
    l = len(PhaseError) * len(phi_mode)
    basic.printProgressBar(0, l, length = 50)
    
    # running index
    index = 0
    
    for indexPhimode, PhaseMode in enumerate(phi_mode):
        
        # Solution for this Phase
        SolutionPhase = []
        

        for indexPhase, PhaseDelay in enumerate(PhaseError): 
           
            # calculate every 100th network stability
            calc_networkstab = True if (indexPhase%Nqy_skip == 0) else False
            
            # calculate in-phase
            SolutionPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                          K_LF = K_LF,
                                          mode=PhaseMode, coupling_function = coupling_function,
                                          Nyquist_Freq = Nqy_freq,
                                          calc_networkstab = calc_networkstab,
                                          feedback_phase = fb_phase,
                                          Export_OLTF=Export_OLTF))
                    
            # Update Progress Bar
            index += 1
            basic.printProgressBar(index, l, length = 50)
            
        # Save Solution
        Solution[str(PhaseMode/np.pi)] = {k: [dic[k] for dic in SolutionPhase] for k in SolutionPhase[0]}
             
    # Some Infos
    elapsed_time = time.time() - start_time
    print("\033[0m" + "Calculations \033[32;1m" + "done" + "\033[0m! (Time needed: " + str(np.round(elapsed_time,3)) + "sec.)")
    
    # Generate VCO Frequency by Scaling with Division Factor
    for SolutionPhase in Solution:
        Solution[SolutionPhase]["VCO_Freq_Stable"] = np.asarray(Solution[SolutionPhase]["Network_Freq_Stable"])*Div_N
        Solution[SolutionPhase]["VCO_Freq_UnStable"] = np.asarray(Solution[SolutionPhase]["Network_Freq_UnStable"])*Div_N   
  
    
    # Some Infos
    print("\033[0m" + "Rescalings \033[32;1m" + "done" + "\033[0m!")
    
    #############################################################################    
    
    return Solution

#############################################################################
###              3th/3.5th Generation, Homogen and Non-Linear
#############################################################################   
def Calc_2Cplg_nonlinear(time_start=0, time_end=10e-9, time_points=1e3,
                      VCO_freq0 = 2*np.pi*24.25e9,
                      Div_N=16*32, K_VCO=2*np.pi*757.47e6,
                      K_PD = 1.6, K_LF = 1,
                      fb_phase = 0,
                      Export_OLTF=False, fitmentset=None):
    #############################################################################    
    """
    Simple Network of two coupled PLL System in Laplace Domain with default values from
    3. Gen Prototypes (Version 1.1) with 2nd order RC-Loopfilter
    
    parameters              description
    =====================  =============================================:
    phase_start             starting phase for iteration
    phase_end               ending phase for iteration
    phase_points            number of iterations
    
    VCO_freq0               VCO Closed Loop Freq (Hz) (default: 24.25e9)
    Div_N                   Divider Value (default: 16*32)
    K_VCO                   (not used) VCO Sensitivity (Hz/V) (default: 757e6)
    K_PD                    (not used) PD Sensitivity (V/pi) (default: 1.6)
    K_LF                    (not used) Loop Filter Gain (default: 1, no laplace)
       
    fb_phase               (optionally) Feedback Phase
    """
    #############################################################################
    
    # Some Infos
    print("\033[34;1m" + "Simulation of two coupled Oscillators (Non-Linear)")
    print("\033[34;1m" + "-> Starting Time {}s to {}s using {} points".format(basic.EngNot(time_start), basic.EngNot(time_end), basic.EngNot(time_points)))
    print("\033[34;1m" + "-> Closed Loop Freq. {}Hz and Division by {}\033[0m".format(basic.EngNot(VCO_freq0/(2*np.pi)), Div_N))
  
    start_time = time.time()
    
    #############################################################################
    #           Variables of PLL System and Network Analysis (Default Values)
    #############################################################################
    
    # laplace variable
    s = sympy.Symbol('s', complex=True)
    
    # Quiescent PLL Freq
    Omega0_Div = VCO_freq0/Div_N

    # Time Delay
    Time_Delay = np.linspace(time_start,time_end,int(time_points))
    
    # PD PETF
    coupling_function = "triangle"
    
    # Frequency Axis for Nyquist Plot
    Nqy_freq = np.linspace(1e3, 100e6, int(1e3))
    
    #############################################################################
    #           Nonlinear Function (including derivertive)
    #############################################################################
    
    # transfer function of PD Input after PETF (Xpd) to CrossCoupling output in steady state
    def  G_CPLG(H_PD, H_prime_PD=None, derivative=False):
        
        if fitmentset == None:
            
            # fitment function parameters (from VCO_FreqVolt fit, ISCAS paper)
            slope = 2.19447032e9 #1.94522384e9
            VCO_qp = 20.80e9 #21.2284874e9
            Vprebias = 2.45 #2.55
            AdderGain = 1.0 #1.2 #1            
        
        else:
            slope, VCO_qp, Vprebias, AdderGain = fitmentset
            
        # fitment function parameters (from VCO_FreqVolt fit, ISCAS paper)

        
        # Calulate Network frequency based on normalized PD Output
        if not(derivative):
            
            # calculate PD output Voltage, based on H_PD and Offset of VCO Pre-Biasffset
            A_PD = AdderGain * 2 * 0.8 * H_PD 
            
            # calculate vco frequency with nonlinear function
            f_VCO = 2 * np.pi *( VCO_qp + slope * np.sqrt(Vprebias + A_PD) )
            
            # calulate network frequency
            f_NET = f_VCO * 1/Div_N
            
            return f_NET

            
        # Calulate stability based on normalized PD Output (Both H_PD inputs are needed)
        else:
            
            # calculate PD output Voltage, based on H_PD and Offset of VCO Pre-Biasffset
            A_PD = AdderGain * 2 * 0.8 * H_PD
            
            # use in steady state this voltage to calcululate vco sensitivity (kvco)
            K_VCO = (2*np.pi*slope) / (2*np.sqrt(Vprebias + A_PD))
 
            # calulate Coupling Gain based on non-linear frequency sensitivity
            CPLG_Gain = H_prime_PD * K_VCO * 1/Div_N
            
            return CPLG_Gain

    #############################################################################
    #           Calculate PLL Solution and Network Stability
    #############################################################################
    
    # Mean Phase for Time Delay
    PhaseError = Omega0_Div * Time_Delay
    
    # Correct Phase Delay by shift
    PhaseError =  PhaseError - fb_phase
    
         
    Solution_InPhase = []
    Solution_AntiPhase = []

    # Progressbar
    l = len(PhaseError)
    basic.printProgressBar(0, l, length = 50)
    
    for index, PhaseDelay in enumerate(PhaseError): 
       
        # calculate every 100th network stability
        calc_networkstab = True if (index%150 == 0) else False
        
        # calculate in-phase
        Solution_InPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                              mode=0, nonlinear=True, K_LF=K_LF,
                                              coupling_function = coupling_function,
                                              Nyquist_Freq = Nqy_freq,
                                              calc_networkstab=calc_networkstab,
                                              Export_OLTF=Export_OLTF))
        
        # calculate anti-phase    
        Solution_AntiPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                                mode=np.pi, nonlinear=True, K_LF=K_LF,
                                                coupling_function = coupling_function,
                                                Nyquist_Freq = Nqy_freq,
                                                calc_networkstab=calc_networkstab,
                                                Export_OLTF=Export_OLTF))
    
 
       # Update Progress Bar
        basic.printProgressBar(index + 1, l, length = 50) 
        
    # Convert List with Dicts to Dicts with Lists
    Solution_InPhase = {k: [dic[k] for dic in Solution_InPhase] for k in Solution_InPhase[0]}
    Solution_AntiPhase = {k: [dic[k] for dic in Solution_AntiPhase] for k in Solution_AntiPhase[0]}
    
    # Some Infos
    elapsed_time = time.time() - start_time
    print("\033[0m" + "Calculations \033[32;1m" + "done" + "\033[0m! (Time needed: " + str(np.round(elapsed_time,3)) + "sec.)")
    
    # Generate VCO Frequency by Scaling with Division Factor
    Solution_InPhase["VCO_Freq_Stable"] = np.asarray(Solution_InPhase["Network_Freq_Stable"])*Div_N
    Solution_InPhase["VCO_Freq_UnStable"] = np.asarray(Solution_InPhase["Network_Freq_UnStable"])*Div_N   
    Solution_AntiPhase["VCO_Freq_Stable"] = np.asarray(Solution_AntiPhase["Network_Freq_Stable"])*Div_N
    Solution_AntiPhase["VCO_Freq_UnStable"] = np.asarray(Solution_AntiPhase["Network_Freq_UnStable"])*Div_N   
    
    # Some Infos
    print("\033[0m" + "Post-Processing \033[32;1m" + "done" + "\033[0m!")
    
    #############################################################################    
    
    return {"m_0": Solution_InPhase,  "m_pi": Solution_AntiPhase}

def Calc_nonlinear(time_start=0, time_end=10e-9, time_points=1e3,
                      VCO_freq0 = 2*np.pi*24.25e9, s = sympy.Symbol('s', complex=True),
                      Div_N=16*32, K_VCO=2*np.pi*757.47e6,
                      K_PD = 1.6, K_LF = 1, Nqy_skip = 10, 
                      fb_phase = 0, phi_mode = [0, np.pi],
                      fitmentset=None, Export_OLTF=False):
    #############################################################################    
    """
    Simple Network of two coupled PLL System in Laplace Domain with default values from
    3. Gen Prototypes (Version 1.1) with 2nd order RC-Loopfilter
    
    parameters              description
    =====================  =============================================:
    phase_start             starting phase for iteration
    phase_end               ending phase for iteration
    phase_points            number of iterations
    
    VCO_freq0               VCO Closed Loop Freq (Hz) (default: 24.25e9)
    Div_N                   Divider Value (default: 16*32)
    K_VCO                   (not used) VCO Sensitivity (Hz/V) (default: 757e6)
    K_PD                    (not used) PD Sensitivity (V/pi) (default: 1.6)
    K_LF                    (not used) Loop Filter Gain (default: 1, no laplace)
       
    fb_phase               (optionally) Feedback Phase
    """
    #############################################################################
    
    # Some Infos
    print("\033[34;1m" + "Simulation of two coupled Oscillators (Non-Linear)")
    print("\033[34;1m" + "-> Starting Time {}s to {}s using {} points".format(basic.EngNot(time_start), basic.EngNot(time_end), basic.EngNot(time_points)))
    print("\033[34;1m" + "-> Closed Loop Freq. {}Hz and Division by {}\033[0m".format(basic.EngNot(VCO_freq0/(2*np.pi)), Div_N))
  
    start_time = time.time()
    
    #############################################################################
    #           Variables of PLL System and Network Analysis (Default Values)
    #############################################################################
    
    # Quiescent PLL Freq
    Omega0_Div = VCO_freq0/Div_N

    # Time Delay
    Time_Delay = np.linspace(time_start,time_end,int(time_points))
    
    # PD PETF
    coupling_function = "triangle"
    
    # Frequency Axis for Nyquist Plot
    Nqy_freq = np.linspace(1e3, 100e6, int(1e3))
    
    #############################################################################
    #           Nonlinear Function (including derivertive)
    #############################################################################
    
    # transfer function of PD Input after PETF (Xpd) to CrossCoupling output in steady state
    def  G_CPLG(H_PD, H_prime_PD=None, derivative=False):
        
        if fitmentset == None:
            
            # fitment function parameters (from VCO_FreqVolt fit, ISCAS paper)
            slope = 1.94522384e9 #2.19447032e9 #1.94522384e9
            VCO_qp = 20.80e9 #21.2284874e9
            Vprebias = 1.95 #2.55 #2.45 #2.55
            AdderGain = 1.1 #1.2 #1            
        
        else:
            slope, VCO_qp, Vprebias, AdderGain = fitmentset
        
        # Calulate Network frequency based on normalized PD Output
        if not(derivative):
            
            # calculate PD output Voltage, based on H_PD and Offset of VCO Pre-Biasffset
            A_PD = AdderGain * 2 * 0.8 * H_PD 
            
            # calculate vco frequency with nonlinear function
            f_VCO = 2 * np.pi *( VCO_qp + slope * np.sqrt(Vprebias + A_PD) )
            
            # calulate network frequency
            f_NET = f_VCO * 1/Div_N
            
            return f_NET

            
        # Calulate stability based on normalized PD Output (Both H_PD inputs are needed)
        else:
            
            # calculate PD output Voltage, based on H_PD and Offset of VCO Pre-Biasffset
            A_PD = 2 * 0.8 * H_PD
            
            # use in steady state this voltage to calcululate vco sensitivity (kvco)
            K_VCO = (2*np.pi*slope) / (2*np.sqrt(Vprebias + A_PD))
 
            # calulate Coupling Gain based on non-linear frequency sensitivity
            CPLG_Gain = H_prime_PD * K_VCO * 1/Div_N
            
            return CPLG_Gain

    #############################################################################
    #           Calculate PLL Solution and Network Stability
    #############################################################################
    
    # Mean Phase for Time Delay
    PhaseError = Omega0_Div * Time_Delay
    
    # Correct Phase Delay by shift
    PhaseError =  PhaseError - fb_phase 
         
    Solution = {}
    
    # Progressbar
    l = len(PhaseError)
    basic.printProgressBar(0, l, length = 50)
 
    # Progressbar
    l = len(PhaseError) * len(phi_mode)
    basic.printProgressBar(0, l, length = 50)
    
    # running index
    index = 0
    
    for indexPhimode, PhaseMode in enumerate(phi_mode):
        
        # Solution for this Phase
        SolutionPhase = []
        

        for indexPhase, PhaseDelay in enumerate(PhaseError): 
           
            # calculate every 100th network stability
            calc_networkstab = True if (indexPhase%Nqy_skip == 0) else False
            
            # calculate in-phase
            SolutionPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                          mode=PhaseMode, nonlinear=True,
                                          coupling_function = coupling_function,
                                          Nyquist_Freq = Nqy_freq,K_LF=K_LF,
                                          calc_networkstab = calc_networkstab,
                                          feedback_phase = fb_phase,
                                          Export_OLTF=Export_OLTF))
                    
            # Update Progress Bar
            index += 1
            basic.printProgressBar(index, l, length = 50)
            
        # Save Solution
        Solution[str(PhaseMode/np.pi)] = {k: [dic[k] for dic in SolutionPhase] for k in SolutionPhase[0]}
         
     
    # Some Infos
    elapsed_time = time.time() - start_time
    print("\033[0m" + "Calculations \033[32;1m" + "done" + "\033[0m! (Time needed: " + str(np.round(elapsed_time,3)) + "sec.)")
    
    # Generate VCO Frequency by Scaling with Division Factor
    for SolutionPhase in Solution:
        Solution[SolutionPhase]["VCO_Freq_Stable"] = np.asarray(Solution[SolutionPhase]["Network_Freq_Stable"])*Div_N
        Solution[SolutionPhase]["VCO_Freq_UnStable"] = np.asarray(Solution[SolutionPhase]["Network_Freq_UnStable"])*Div_N   
      
    # Some Infos
    print("\033[0m" + "Post-Processing \033[32;1m" + "done" + "\033[0m!")
    
    #############################################################################    
    
    return Solution

#############################################################################
###                 Dynamical Model
#############################################################################
  
# Simulation of coupled oscillators
def Run_Timeseries(Osc_Param, Network_Param, tstop, dt,
                   threshold_steadystate_phase=1e-4,
                   threshold_steadystate_freq=5e-5,
                   PrintOutput=True, EnableCouplingAtDelay=False):
    
    #############################################################################    
    """
    Calculate Coupled Networks using Euler Scheme
    
    parameters              description
    =====================  =============================================:
    tstop                   Stop time or maximale time
    dt                      time step
    
    Osc_Param               array containing oscillator parameters as dict
    
                            Osc_Param[<osc_index>]["omega0"]		# natural frequency
                            Osc_Param[<osc_index>]["theta0"]		# initial phase
                            Osc_Param[<osc_index>]["K_vco"]			# Sensitivity of the VCO
							Osc_Param[<osc_index>]["pd_amplitude"]	# Amplitude of PD output
							Osc_Param[<osc_index>]["pd_fct"]		# PD function: XOR or SIN
							Osc_Param[<osc_index>]["lf_fct"]		# LF function: 1_RC
                            Osc_Param[<osc_index>]["lf_gain"]		# LF gain (only of lf_fct is defined)              
                            Osc_Param[<osc_index>]["lf_tau"]		# LF time constants (only of lf_fct is defined)   
                            
    Network_Param           dictionary containing network parameters

                            Network_Param["cplg_gain"]  = < matrix >
                            Network_Param["cplg_delay"] = < matrix >
                            
                            
    threshold_steadystate_phase         threshold for detection of steady state of the phase
    threshold_steadystate_freq          threshold for detection of steady state of the frequency 
                      
    EnableCouplingAtDelay               disable coupling until delay is reached
                            
    example with 3 nodes in chain


            Omega0_Net = 2 * np.pi * 23.5e9 / 512
            Kpd = 1.6/2
            lf_tau = 220*680e-12
            Kvco = (2*np.pi*757.47e6)/512
            
            delay1 = 20e-9
            delay2 = 20e-9
            
            Osc_Param = []
            
            # oscillator A
            Osc_Param.append( {"omega0": Omega0_Net,
                               "theta0": 0,
							   "pd_amplitude": Kpd,
							   "pd_fct": "SIN",
							   "K_vco": Kvco,
							   "lf_gain":1,
							   "lf_fct": "RC_1",
							   "lf_tau": lf_tau } )
            
            # oscillator B
            Osc_Param.append( {"omega0": Omega0_Net,
                               "theta0": 1.1*np.pi,
							   "pd_amplitude": Kpd,
							   "pd_fct": "SIN",
							   "K_vco": Kvco,
							   "lf_gain":1,
							   "lf_fct": "RC_1",
							   "lf_tau": lf_tau } )
            
            # oscillator C
            Osc_Param.append( {"omega0": Omega0_Net,
                               "theta0": 0.9*np.pi,
							   "pd_amplitude": Kpd,
							   "pd_fct": "SIN",
							   "K_vco": Kvco,
							   "lf_gain":1,
							   "lf_fct": "RC_1",
							   "lf_tau": lf_tau } )
            
            # Network Topology
            Network_Param = {
                            "cplg_gain": np.matrix([[0, 1, 0],
                                                    [1, 0, 1],
                                                    [0, 1, 0]]),
                            
                             "cplg_delay": np.matrix([[0, delay1, 0],
                                                      [delay1, 0, delay2],
                                                      [0, delay2, 0]])
                             }
            
            
            [time, results, info] = sll.Run_Timeseries(Osc_Param,Network_Param, 2e-6, 1e-9) 
    """
    #############################################################################

    # Steady State detector
    def SteadyStateDetection(dataset, start_index, threshold_steadystate):
 
        # Find number of threshold times true in a row
        def FindTrue(arr, threshold = 25):
            count = 0
    
            for i, val in enumerate(arr):
                if val: count += 1
                else: count = 0
                    
                if count == threshold: return i - threshold + 1
                
            return -1
        
        Delta_Norm = (dataset - dataset[-1])
        Delta_Threshold = ( np.abs( np.diff(Delta_Norm[start_index::]) ) <  threshold_steadystate )[:-1]
        
        Index_SteadyState = FindTrue( Delta_Threshold )
        
        # Check if negative
        if Index_SteadyState < 0: Index_SteadyState = -1    
        
        # Return Index
        return Index_SteadyState        

    #############################################################################
    
    # return dictionary
    results = []

    # simulation time series
    t = np.arange(0, tstop, dt)

    # Some Infos
    if PrintOutput:
        print("\033[34;1m" + "Time Series Simulation of coupled Oscillators")
        print("\033[34;1m" + "-> Stop Time {}s using step size {}s ({} points)".format(basic.EngNot(tstop),
                                                                                       basic.EngNot(dt),
                                                                                       basic.EngNot(len(t))) )
    start_time = time.time()

    #############################################################################
    
    # Initialize solution arrays
    for osc_index in range(len(Osc_Param)):
        
        oscdata = {}
        
        # node voltages / values
        oscdata["vadd"]             = np.zeros(len(t))      # adder voltage
        oscdata["vtune"]            = np.zeros(len(t))      # tuning voltage
        oscdata["theta"]             = np.zeros(len(t))     # oscillator instantaneous phase
        oscdata["angular_frequency"] = np.zeros(len(t))     # oscillator frequency
                      

        # init values for time 0
        oscdata["angular_frequency"][0] = Osc_Param[osc_index]["omega0"]
        oscdata["theta"][0] = Osc_Param[osc_index]["theta0"]
            
        results.append(oscdata)       
 
    
 
        # Info - Node Settings
        s_nodes = "\033[34;1m"
        s_nodes += "-> Oscillator {}: natural frequency {}Hz ({}radHz) and initial phase {}rad".format(osc_index+1,
                                                                                                       basic.EngNot(Osc_Param[osc_index]["omega0"]/(2*np.pi)),
                                                                                                       basic.EngNot(Osc_Param[osc_index]["omega0"]),
                                                                                                       basic.EngNot(Osc_Param[osc_index]["theta0"]) )
        s_nodes += "\033[0m"
        

        
        # This Oscillator is coupled to...
        s_coupling ="\033[34;1m" + "   Couplings: "
        
        # Coupled Oscillators to this Index
        for cplg_index in np.transpose( np.nonzero( Network_Param["cplg_gain"][osc_index, ::] ) ):
            s_coupling += "{} <- {}s -> {} ".format(osc_index+1,
                                                  basic.EngNot(Network_Param["cplg_delay"][osc_index, cplg_index[1]]),
                                                  cplg_index[1]+1)
        
        s_coupling += "\033[0m"

        if PrintOutput:
            print( s_nodes )        
            print( s_coupling )

    #############################################################################
        
    # Progressbar
    basic.printProgressBar(0, len(t), length = 50)
        
    # calculus
    for i in range(1, len(t)):
        
        # Calulation
        Perform_Timeseries_Simulation(i, dt, results, Osc_Param, Network_Param, CheckEnableCplg=EnableCouplingAtDelay)
        
        # Update Progress Bar
        basic.printProgressBar(i + 1, len(t), length = 50) 
        
 
    # Some Infos
    elapsed_time = time.time() - start_time
    if PrintOutput:
        print("\033[0m" + "Calculations \033[32;1m" + "done" + "\033[0m! (Time needed: " + str(np.round(elapsed_time,3)) + "sec.)")
        
    #############################################################################
    # Post Process each node
    for osc_index in range(len(Osc_Param)):
        
        # Post Process Data
        results[osc_index]["frequency"] =  results[osc_index]["angular_frequency"]/(2*np.pi)
        results[osc_index]["phase"] =  results[osc_index]["theta"] * 180/np.pi
        
        # Info - Node
        s_nodes = "\033[34;1m"
        s_nodes += "-> Oscillator {}: final frequency {}Hz".format(osc_index+1, basic.EngNot(results[osc_index]["frequency"][-1], sig_figs=5) )
        s_nodes += "\033[0m"        
        
        if PrintOutput:
            print( s_nodes )         
        
    if PrintOutput:
        print("\033[0m" + "Post Processing of Node Parameters \033[32;1m" + "done" + "\033[0m!")        

    #############################################################################

    postprocess_coupling = []

    # Post Process Coupling
    for osc_index in range(len(Osc_Param)):
        
        
        # Iterate each coupled Oscillator
        for cplg_index in np.transpose( np.nonzero( Network_Param["cplg_gain"][osc_index, ::] ) ):  
            
            cplg_result = {}
            cplg_result["cplg"] = [osc_index+1, cplg_index[1]+1]
            cplg_result["success"] = True
            
            # Extract Deltas
            cplg_result["delta_theta"] = results[osc_index]["theta"] - results[cplg_index[1]]["theta"]
            cplg_result["delta_phase"] =  cplg_result["delta_theta"] * 180/np.pi
            cplg_result["delta_frequency"] = results[osc_index]["frequency"] - results[cplg_index[1]]["frequency"]
            
            # index of coupling start
            cplg_start = int(Network_Param["cplg_delay"][osc_index, cplg_index[1]]/dt)
            
            # Find settling time
            cplg_result["settlingtime_theta"] = t[ SteadyStateDetection(cplg_result["delta_theta"], cplg_start, threshold_steadystate_phase) ]
            cplg_result["settlingtime_frequency"] = t[ SteadyStateDetection(cplg_result["delta_frequency"], cplg_start, threshold_steadystate_freq) ]
 
            cplg_result["settlingtime"] = np.mean([cplg_result["settlingtime_theta"], cplg_result["settlingtime_frequency"]])
            
            # Check if settling time is smaller than time delay
            if( cplg_result["settlingtime"] >= t[-2]):
                print ("\033[1;31mWARNING - Coupling Path: {}<->{} not settled!\033[0m".format(osc_index+1, cplg_index[1]+1))
                cplg_result["success"] = False

            # Check if settling time is before coupling stated
            if( cplg_result["settlingtime"] <= t[ cplg_start ]):
                print ("\033[1;31mWARNING - Coupling Path: {}<->{} settled before activated!\033[0m".format(osc_index+1, cplg_index[1]+1))
                cplg_result["success"] = False
                
            # initial phase difference at active coupling
            cplg_result["initial_delta_theta"] = cplg_result["delta_theta"][ cplg_start ]
            
            postprocess_coupling.append(cplg_result)

            # Info - Coupling
            s_coupling = "\033[34;1m"
            s_coupling += "-> Coupling Path: {}<->{}: final frequency: {}Hz (delta: {}Hz) final phase difference {}rad".format(osc_index+1,
                                                                                                                               cplg_index[1]+1,
                                                                                                                               basic.EngNot(results[osc_index]["frequency"][-1], sig_figs=5),
                                                                                                                               basic.EngNot(cplg_result["delta_frequency"][-1], sig_figs=5),
                                                                                                                               basic.EngNot(cplg_result["delta_theta"][-1], sig_figs=5))
            s_coupling += "\033[0m"        
            
            if PrintOutput:
                print( s_coupling )  
                
    if PrintOutput:
        print("\033[0m" + "Post Processing of Coupling Parameters \033[32;1m" + "done" + "\033[0m!")

    #############################################################################
              
    return [t, results, postprocess_coupling]
    
