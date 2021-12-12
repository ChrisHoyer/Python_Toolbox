#############################################################################
#   Different Scripts to Model State-Locked-Loops
#
#   - Net_2Mutually: Simple Network in Laplace with two bidirectional coupled nodes to Calc Global Freq
#
#   Prototype Models
#   - Calc_2Cplg_linear: Calculates important Values for Steady State Analysis
#   - Calc_2Cplg_nonlinear: Calculates important Values for Steady State Analysis (nonlinear)
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

#############################################################################
###                 Linear Network with 2 Mutually Coupled Oscillators
#############################################################################
def Net_2Mutually(phase_delay, omega0_div, G_CPLG, variable, 
                  mode = 0, nonlinear = False,
                  feedback_phase = 0,
                  coupling_function= "triangle",
                  Nyquist_Freq = np.linspace(0,1e6,int(1e5)),
                  G_CPLG_LF = None,
                  Nyquist_Neglect = 10,
                  calc_networkstab = True,
                  nyquist_tolerance = 1e-3,
                  stab_tolerance = 1e-14,
                  debug_return = True, freqinHz = True):
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
    G_CPLG_LF               (optionally) G_CPLG but for Network  with LF (no steady state)
    nyquist_tolerance       (optionally) tolerance of nyquist solving
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
        
    # Use G_CPLG also for Network Stability
    if G_CPLG_LF == None:
        G_CPLG_LF = G_CPLG
    
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
 
# =========================================================================    
       
    else:
        
        # Transfer Function of Open Loop with H_prime
        G_CPLG_TF = G_CPLG * H_prime_PD

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
#           Calculate Network Stability
#############################################################################  
    
    # Export no Values:
    return_var["Nyquist_Point"] =  float("NaN")  
    
    # return Nyquist Calculation
    if debug_return:
        return_var["Nyquist_Solution"] = float("NaN")
        return_var["Nyquist_Freq"] = float("NaN")

            
    # Calculate Network Stability
    if calc_networkstab and not(nonlinear):
        
        # Laplace Delay block for calculacted Delay
        exp_delay = sympy.exp(-time_delay*variable) 
        
        # Transfer Function of Open Loop with LF with H_prime
        G_CPLG_TF_LF = G_CPLG_LF * H_prime_PD
    
        # Closed Loop Transfer Function of a Single PLL with LF
        G_CPLG_CL_LF = G_CPLG_TF_LF / (1 + G_CPLG_TF_LF)
        G_CPLG_CL_LF = sympy.simplify(G_CPLG_CL_LF)
        
        # Open loop network transfer function of coupled PLLs
        G_NET = G_CPLG_CL_LF * exp_delay * G_CPLG_CL_LF * exp_delay
        
        # Calculate Nyquist Frequency Span
        Nyquist_Omega = 2j * np.pi * Nyquist_Freq
            
        # Nyquist Plot with Complex return
        Nyquist_Calc = ctrl.Extract_Sympy_1Var(G_NET, Nyquist_Omega,
                                               variable=variable)
        
        # return Nyquist Calculation
        if debug_return:
            return_var["Nyquist_Solution"] = Nyquist_Calc
            return_var["Nyquist_Freq"] = Nyquist_Omega*scaletoHz
                  
        # Only one Value? Or no Max Value
        try:
            
            # Find Entry Point in Quadrant 1, where Imag and Real Part is positive
            index_NyquistPoint_Q1 = np.where((np.real(Nyquist_Calc) > 0) & (np.imag(Nyquist_Calc) > 0))
            index_NyquistPoint_Q1 = np.min(index_NyquistPoint_Q1)
        
            # Find Values, where the imag part is close to zero
            index_NyquistPoint = np.argmin(np.abs(np.imag(Nyquist_Calc[index_NyquistPoint_Q1:-1])))
        
            # Save Nyquist Point and Distance from Real Part to 1
            return_var["Nyquist_Point"] = Nyquist_Calc[index_NyquistPoint_Q1+index_NyquistPoint]

        except:
                
            #print("An exception occurred")
            pass
            
    return return_var


#############################################################################
###                 Different Prototype Routines
#############################################################################

                                  
#############################################################################
###              3th/3.5th Generation, Homogen and Linear
#############################################################################  
def Calc_2Cplg_linear(time_start=0, time_end=10e-9, time_points=1e3,
                      VCO_freq0 = 2*np.pi*24.25e9,
                      Div_N=16*32, K_VCO=2*np.pi*757.47e6,
                      K_PD = 1.6, K_LF = 1, fb_phase = 0,  
                      s = sympy.Symbol('s'), Nqy_freq = np.linspace(1e3, 100e6, int(1e3)),
                      Nqy_skip = 10):
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
    G_CPLG_LF = K_PD * K_LF  * K_VCO/s * 1/Div_N
    
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
                                               G_CPLG_LF = G_CPLG_LF,
                                               mode=0, coupling_function = coupling_function,
                                               Nyquist_Freq = Nqy_freq,
                                               calc_networkstab = calc_networkstab,
                                               feedback_phase = fb_phase))
        
        # calculate anti-phase    
        Solution_AntiPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                                 G_CPLG_LF = G_CPLG_LF,
                                                 mode=np.pi, coupling_function = coupling_function,
                                                 Nyquist_Freq = Nqy_freq,
                                                 calc_networkstab = calc_networkstab,
                                                 feedback_phase = fb_phase))
        
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

#############################################################################
###              3th/3.5th Generation, Homogen and Non-Linear
#############################################################################   
def Calc_2Cplg_nonlinear(time_start=0, time_end=10e-9, time_points=1e3,
                      VCO_freq0 = 2*np.pi*24.25e9,
                      Div_N=16*32, K_VCO=2*np.pi*757.47e6,
                      K_PD = 1.6, K_LF = 1,
                      fb_phase = 0 ):
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
    s = sympy.Symbol('s')
    
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
        
        # fitment function parameters (from VCO_FreqVolt fit, ISCAS paper)
        slope = 2.19447032e9 #1.94522384e9
        VCO_qp = 20.80e9 #21.2284874e9
        Vprebias = 2.45 #2.55
        
        # Calulate Network frequency based on normalized PD Output
        if not(derivative):
            
            # calculate PD output Voltage, based on H_PD and Offset of VCO Pre-Biasffset
            A_PD = 2 * 0.8 * H_PD 
            
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
                                              mode=0, nonlinear=True,
                                              coupling_function = coupling_function,
                                              Nyquist_Freq = Nqy_freq,
                                              calc_networkstab=calc_networkstab))
        
        # calculate anti-phase    
        Solution_AntiPhase.append(Net_2Mutually(PhaseDelay, Omega0_Div, G_CPLG, s,
                                                mode=np.pi, nonlinear=True,
                                                coupling_function = coupling_function,
                                                Nyquist_Freq = Nqy_freq,
                                                calc_networkstab=calc_networkstab))
    
 
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


#############################################################################
###              3th/3.5th Generation, Homogen and Non-Linear
#############################################################################   
