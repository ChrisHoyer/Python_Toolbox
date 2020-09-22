#############################################################################
#   Different Scripts to Model State-Locked-Loops
#
#   - Fnet_Bi2N_Stab: Simple Network in Laplace with two bidirectional coupled nodes to Calc Global Freq
#
#   Prototype Models
#   - Calc_3thGen_V1_1_homogen
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Stand: 31.08.2020
#############################################################################
import controltheory_toolbox as ctrl

import sympy
import time
import numpy as np
import scipy as sp
from scipy import signal

#############################################################################
###                 Linear Fnet_Bi2N incl. Stability
#############################################################################
def Fnet_Bi2N_Stab(phase_delay, omega0_div, G_CPLG, variable, 
                   mode = 0, coupling_function= "triangle",
                   Nyquist_Freq = np.linspace(0,1e6,int(1e5)),G_CPLG_LF = None,
                   Nyquist_Neglect = 10, calc_networkstab = True,
                   nyquist_tolerance = 1e-3, stab_tolerance = 1e-14,
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
    variable                laplace variable
    
    mode                    (optionally) mode locking (0 or pi)
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
    phase_error = mode - phase_delay

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
    
    # Locked State Gain (Steady State)
    G_LS = variable * G_CPLG
    G_LS = float(sympy.limit(G_LS, variable, 0))
    
    # calculate global frequency
    OmegaNet = G_LS * H_PD + omega0_div
    return_var["Network_Freq"] = OmegaNet*scaletoHz
       
    # calculate time delay
    time_delay = phase_delay/OmegaNet
    return_var["Time_Delay"] = time_delay
    
#############################################################################   
#           Calculate PLL Stability
#############################################################################  
    
    # Transfer Function of Open Loop with H_prime
    G_CPLG_TF = G_CPLG * H_prime_PD
    
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
    return_var["Nyquist_Gain"] =  float("NaN")  
    
    # return Nyquist Calculation
    if debug_return:
        return_var["Nyquist_Solution"] = float("NaN")
        return_var["Nyquist_Freq"] = float("NaN")

            
    # Calculate Network Stability
    if calc_networkstab:
        
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
            
            # Remove the first 100 points (near Nyquist Point) and get abs value
            Imag_Abs = np.abs(np.imag(Nyquist_Calc))[Nyquist_Neglect:]
        
            # Find Point index where Imaginary Part is Zero (< Tolerance)
            Imag_Zero_index = np.where(Imag_Abs < nyquist_tolerance)[0]
    
            # get corresponding real part
            RealParts = np.real(Nyquist_Calc)[Imag_Zero_index]
            MaxReal_index = np.argmax(RealParts)
        
            # get Peak index
            NyquistPoint_index = Imag_Zero_index[MaxReal_index] 
        
            # Get Max Peak of realpart, when imaginary part is near zero       
            return_var["Nyquist_Gain"] =  np.real(Nyquist_Calc[NyquistPoint_index])
            
        except:
                
            #print("An exception occurred")
            pass
            
    return return_var

#############################################################################
###                 Nonlinear Fnet_Bi2N incl. Stability
#############################################################################
def Fnet_Bi2N_Stab_NonLinear(phase_delay, omega0_div, G_CPLG,  variable,
                   mode = 0, coupling_function= "triangle",
                   Nyquist_Freq = np.linspace(0,1e6,int(1e5)),G_CPLG_LF = None,
                   Nyquist_Neglect = 10, calc_networkstab = False,
                   nyquist_tolerance = 1e-3, stab_tolerance = 1e-14,
                   debug_return = True, freqinHz = True):
#############################################################################    
    """
   Simple Network in Laplace with two bidirectional coupled nodes to calculate
   the global network frequency and stability vs the phase delay 
   Note, that this is a differential equation. Where the init value has to be set.
    
    parameters              description
    =====================  =============================================:
    phase_delay             current phase delay
    omega0_div              (not used) divided quiescent PLL frequency
    G_CPLG                  Nonlineare Open Loop Transfer Function in steady state (K_LF=1)
    variable                laplace variable
    
    mode                    (optionally) mode locking (0 or pi)
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
        
    
    # return dictionary
    return_var = {}
    return_var["Phase_Delay"] = -phase_delay/np.pi
    return_var["Net_Calc"] = calc_networkstab
        
    # Calculate Phase Error tilde substitude
    phase_error = mode - phase_delay

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
    
    # calculate global frequency
    OmegaNet = G_CPLG(H_PD)
    return_var["Network_Freq"] = OmegaNet*scaletoHz
       
    # calculate time delay
    time_delay = phase_delay/OmegaNet
    return_var["Time_Delay"] = time_delay
    
#############################################################################   
#           Calculate PLL Stability
#############################################################################  
    
    # Transfer Function of Open Loop (derivative) with H_prime and H
    G_CPLG_TF = G_CPLG(H_PD, X_prime_PD=H_prime_PD, derivative=True) * 1/variable
    
    # Closed Loop Transfer Function of a Single PLL
    G_CPLG_CL = G_CPLG_TF / (1 + G_CPLG_TF)
    G_CPLG_CL = sympy.simplify(G_CPLG_CL)
    
    # Generate Nominator and Denmominator for Analyzes of Poles/Zeros
    G_CPLG_CL_nom, G_CPLG_CL_denom = sympy.fraction(G_CPLG_CL)   

    # export function for Denom and Nom
    G_CPLG_CL_nom = sympy.lambdify([variable], G_CPLG_CL_nom)
    G_CPLG_CL_denom = sympy.lambdify([variable], G_CPLG_CL_denom)
        
    # initial value for solving for stability
    G_LS = G_CPLG(0)
    
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
 
           
    return return_var

#############################################################################
###                 Different Prototype Routines
#############################################################################

                                  
#############################################################################
###                 3th Generation, Homogen and Linear PCB V1.1
#############################################################################
    
def Calc_3thGen_V1_1_homogen(time_start=0, time_end=10e-9, time_points=1e3,
                              VCO_freq0 = 2*np.pi*24.25e9,
                              Div_N=16*32, K_VCO=2*np.pi*757e6,
                              K_PD = 1.6, K_LF = 1,
                              phaseshift=0,
                              laplace_var = ""):
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
    K_VCO                   VCO Sensitivity (Hz/V) (default: 757e6)
    K_PD                    PD Sensitivity (V/pi) (default: 1.6)
    K_LF                    Loop Filter Gain (default: 1, no laplace)
    
    phaseshift              (optionally) small phase correction shift (times pi)
    laplace_var             (optionally) laplace variable (default: "s")
    """
    #############################################################################
    
    # Some Infos
    print("Starting Calculations!")
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
    #           Calculate PLL Solution and Network Stability
    #############################################################################
    
    # transfer function of PD Input to CrossCoupling output
    G_CPLG = K_PD * K_VCO/s * 1/Div_N
    G_CPLG_LF = K_PD * K_LF  * K_VCO/s * 1/Div_N
    
    # Mean Phase for Time Delay
    PhaseError = Omega0_Div * Time_Delay
    
    # Correct Phase Delay by shift
    PhaseError =  PhaseError - phaseshift*np.pi
         
    Solution_InPhase = []
    Solution_AntiPhase = []
    
    for index, PhaseDelay in enumerate(PhaseError): 
       
        # calculate every 100th network stability
        calc_networkstab = True if (index%150 == 0) else False
        
        # calculate in-phase
        Solution_InPhase.append(Fnet_Bi2N_Stab(PhaseDelay, Omega0_Div, G_CPLG, s,
                                               G_CPLG_LF = G_CPLG_LF,
                                               mode=0, coupling_function = coupling_function,
                                               Nyquist_Freq = Nqy_freq,
                                               calc_networkstab=calc_networkstab))
        
        # calculate anti-phase    
        Solution_AntiPhase.append(Fnet_Bi2N_Stab(PhaseDelay, Omega0_Div, G_CPLG, s,
                                                 G_CPLG_LF = G_CPLG_LF,
                                                 mode=np.pi, coupling_function = coupling_function,
                                                 Nyquist_Freq = Nqy_freq,
                                                 calc_networkstab=calc_networkstab))
    
         
    # Convert List with Dicts to Dicts with Lists
    Solution_InPhase = {k: [dic[k] for dic in Solution_InPhase] for k in Solution_InPhase[0]}
    Solution_AntiPhase = {k: [dic[k] for dic in Solution_AntiPhase] for k in Solution_AntiPhase[0]}
    
    # Some Infos
    elapsed_time = time.time() - start_time
    print("Calculations done! (Time needed: " + str(np.round(elapsed_time,3)) + "sec.)")
    
    # Generate VCO Frequency by Scaling with Division Factor
    Solution_InPhase["VCO_Freq_Stable"] = np.asarray(Solution_InPhase["Network_Freq_Stable"])*Div_N
    Solution_InPhase["VCO_Freq_UnStable"] = np.asarray(Solution_InPhase["Network_Freq_UnStable"])*Div_N   
    Solution_AntiPhase["VCO_Freq_Stable"] = np.asarray(Solution_AntiPhase["Network_Freq_Stable"])*Div_N
    Solution_AntiPhase["VCO_Freq_UnStable"] = np.asarray(Solution_AntiPhase["Network_Freq_UnStable"])*Div_N   
    
    # Some Infos
    print("Rescalings done!")
    
    #############################################################################    
    
    return {"m_0": Solution_InPhase,  "m_pi": Solution_AntiPhase}

#############################################################################
###                 3th Generation, Homogen and Non-Linear PCB V1.1
#############################################################################
    
def Calc_3thGen_V1_1_homogen_nonlinear(time_start=0, time_end=10e-9, time_points=1e3,
                                       Div_N=16*32, phaseshift=0):
    #############################################################################    
    """
    Simple Network of two coupled PLL System in Laplace Domain with default values from
    3. Gen Prototypes (Version 1.1) with 2nd order RC-Loopfilter
    
    parameters              description
    =====================  =============================================:
    phase_start             starting phase for iteration
    phase_end               ending phase for iteration
    phase_points            number of iterations
    
    phaseshift              (optionally) small phase correction shift (times pi)
    laplace_var             (optionally) laplace variable (default: "s")
    """
    #############################################################################
    
    # Some Infos
    print("Starting Calculations!")
    start_time = time.time()
    
    #############################################################################
    #           Variables of PLL System and Network Analysis (Default Values)
    #############################################################################
    
    # laplace variable
    s = sympy.Symbol('s')
        
    # Loop Filter Transfer Function
    K_LF = 1
     
    # Quiescent VCO Frequency
    PLL_freq0 = 2*np.pi*20.7923e9
                         
    # Quiescent PLL Network Freq
    Omega0_Div = PLL_freq0/Div_N

    # Time Delay
    Time_Delay = np.linspace(time_start,time_end,int(time_points))
    
    # PD PETF
    coupling_function = "triangle"
    
     # Frequency Axis for Nyquist Plot
    Nqy_freq = np.linspace(1e3, 100e6, int(1e3))
    
    #############################################################################
    #           Calculate PLL Solution and Network Stability
    #############################################################################
    
    # transfer function of PD Input after PETF (Xpd) to CrossCoupling output in steady state
    def  G_CPLG(X_PD, X_prime_PD=None, derivative=False):
        
        # fitment function parameters
        slope = 2.1974e9
        VCO_qp = 20.7923e9
        Vprebias = 2.5
        
        # Calulate Network frequency based on normalized PD Output
        if not(derivative):
            
            # calculate PD output Voltage, based on H_PD and Offset of VCO Pre-Biasffset
            A_PD = 2 * 0.8 * X_PD 
            
            # calcululate vco frequency with nonlinear function
            f_VCO = 2*np.pi (VCO_qp + slope * np.sqrt(Vprebias + A_PD))
            
            # calulate network frequency
            f_NET = f_VCO * 1/Div_N
            
            return f_NET

            
        # Calulate stability based on normalized PD Output
        else:
            
            # calculate PD output Voltage, based on H_PD and Offset of VCO Pre-Biasffset
            A_PD = 2 * 0.8 * X_P
            
            # use in steady state this voltage to calcululate vco sensitivity
            # multiplied by inner function
            K_VCO = X_prime_PD * (2*np.pi*slope) / (2*np.sqrt(Vprebias + A_PD))
 
            # calulate network frequency sensitivity
            f_NET = K_VCO* 1/Div_N
            
            return f_NET

    #############################################################################
    # Mean Phase for Time Delay
    PhaseError = Omega0_Div * Time_Delay
    
    # Correct Phase Delay by shift
    PhaseError =  PhaseError - phaseshift*np.pi
         
    Solution_InPhase = []
    Solution_AntiPhase = []
    
    for index, PhaseDelay in enumerate(PhaseError): 
       
        # calculate every 100th network stability
        calc_networkstab = True if (index%150 == 0) else False
        
        # calculate in-phase
        Solution_InPhase.append(Fnet_Bi2N_Stab_NonLinear(PhaseDelay, Omega0_Div, G_CPLG, s,
                                               mode=0, coupling_function = coupling_function,
                                               Nyquist_Freq = Nqy_freq,
                                               calc_networkstab=calc_networkstab))
        
        # calculate anti-phase    
        Solution_AntiPhase.append(Fnet_Bi2N_Stab_NonLinear(PhaseDelay, Omega0_Div, G_CPLG, s,
                                                 mode=np.pi, coupling_function = coupling_function,
                                                 Nyquist_Freq = Nqy_freq,
                                                 calc_networkstab=calc_networkstab))
    
         
    # Convert List with Dicts to Dicts with Lists
    Solution_InPhase = {k: [dic[k] for dic in Solution_InPhase] for k in Solution_InPhase[0]}
    Solution_AntiPhase = {k: [dic[k] for dic in Solution_AntiPhase] for k in Solution_AntiPhase[0]}
    
    # Some Infos
    elapsed_time = time.time() - start_time
    print("Calculations done! (Time needed: " + str(np.round(elapsed_time,3)) + "sec.)")
    
    # Generate VCO Frequency by Scaling with Division Factor
    Solution_InPhase["VCO_Freq_Stable"] = np.asarray(Solution_InPhase["Network_Freq_Stable"])*Div_N
    Solution_InPhase["VCO_Freq_UnStable"] = np.asarray(Solution_InPhase["Network_Freq_UnStable"])*Div_N   
    Solution_AntiPhase["VCO_Freq_Stable"] = np.asarray(Solution_AntiPhase["Network_Freq_Stable"])*Div_N
    Solution_AntiPhase["VCO_Freq_UnStable"] = np.asarray(Solution_AntiPhase["Network_Freq_UnStable"])*Div_N   
    
    # Some Infos
    print("Rescalings done!")
    
    #############################################################################    
    
    return {"m_0": Solution_InPhase,  "m_pi": Solution_AntiPhase}
