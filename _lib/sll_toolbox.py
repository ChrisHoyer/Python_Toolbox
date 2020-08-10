#############################################################################
#   Different Scripts to Model State-Locked-Loops
#
#   - Fnet_Bi2N: Simple Network in Laplace with two bidirectional coupled nodes to Calc Global Freq
#   - Fnet_Bi2N_Stab: Network of Two Coupled PLLs
#
#   Prototype Models
#   - Model_3thGen_V1_1_homogen
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Stand: 10.01.2020
#############################################################################
import controltheory_toolbox as ctrl
import basic_toolbox as basic


import sympy
import numpy as np
import scipy as sp
from scipy import signal
import matplotlib.pyplot as plt

#############################################################################
###                 Fnet_Bi2N 
#############################################################################
def Fnet_Bi2N(last_state, delay, omega0_div, GTF_FF, GTF_FB, variable, 
              INV = False, phaseshift=-0.5,
              calc_stability = True, stab_tolerance = 1e-14,
              coupling_function= "triangle"):
#############################################################################    
    """
   Simple Network in Laplace with two bidirectional coupled nodes to calculate
   the global network frequency vs the phase delay. Stability for only one note
   is calculated. Both nodes in the network are identical.
    
    parameters              description
    =====================  =============================================:
    last_state              last state from (as tuple!)
    delay                   current delay value
    freq0_div               divided center frequency
    GTF_FF                  Feed Forward Transfer Function (with laplace variable)
    GTF_FB                  Feed Back Transfer Function (with laplace variable)
    variable                laplace variable
    INV                     (optionally) phase shift by pi
    phaseshift              (optionally) phase offset related (in pi)
    calc_stability          (optionally) calculate stability of the plots
    stab_tolerance          (optionally) tolerance of stability solver
    coupling_function       (optionally) PD coupling function (XOR=sawtooth, Mixer=cos)

    return type
       next_state and time delay as tuple [x0, x1, tau0, tau1]
    
    Example:
        
        import sll_toolbox as sll
        

    """
#############################################################################   
       

    # delay is respresented as phase in radian (0... 2pi)
    tau_phase0 = delay
    tau_phase1 = delay
    
    # Calculate phase difference based on delay in radian
    phase_difference0 = last_state[0] - tau_phase0 - last_state[1]
    phase_difference1 = last_state[1] - tau_phase1 - last_state[0]
 
    # Add Invertion and Phaseshifts (Represents varphi_const), which is equal
    phase_difference0 = phase_difference0  - (INV*np.pi) - phaseshift*np.pi
    phase_difference1 = phase_difference1  - (INV*np.pi) - phaseshift*np.pi

    # coupling function, XOR
    if coupling_function == "triangle":
        
        # Calculate Sawthooth function        
        delta_phase0 = 0.5*sp.signal.sawtooth(phase_difference0, width=0.5)
        delta_phase1 = 0.5*sp.signal.sawtooth(phase_difference1, width=0.5)
        
        # derivertive of the slope
        PD_Derivertive = 1/np.pi*sp.signal.square(phase_difference0, duty=0.5)

    # coupling function, MIXER
    elif coupling_function == "cos":
        
        # Calculate Sawthooth function        
        delta_phase0 = -1 * np.cos(phase_difference0)
        delta_phase1 = -1 * np.cos(phase_difference1)
        
        # derivertive of the slope
        PD_Derivertive = np.sin(phase_difference1)
        
    # coupling function, MIXER
    elif coupling_function == "sawtooth":
        
        # Calculate Sawthooth function        
        delta_phase0 = sp.signal.sawtooth(phase_difference0)
        delta_phase1 = sp.signal.sawtooth(phase_difference1)
        
        # derivertive of the slope
        PD_Derivertive = sp.signal.square(phase_difference0, duty=1)  
        
#############################################################################  
    
    # Locked State Gain (Steady State)
    G_LS = variable * GTF_FF * GTF_FB
    G_LS = float(sympy.limit(G_LS, variable, 0))
    
    # calculate global frequency
    omega0 = G_LS * delta_phase0 + omega0_div
    omega1 = G_LS * delta_phase1 + omega0_div
    
    # calculate time delay
    tau0 = delay/omega0
    tau1 = delay/omega1

#############################################################################  
    
    # calculate stability of one single node
    if calc_stability:
        
        # derivertive to PD slope
        GTF_FF_PD = PD_Derivertive * GTF_FF * GTF_FB
                
        # Generate transfer function dependcy on phase difference for single PLL
        G_CPLG = GTF_FF_PD / (1 + GTF_FF_PD)
        
        # Generate Function (Denom, Nom) for coupling, that has to be solved
        G_CPLG = sympy.simplify(G_CPLG)
        G_CPLG_nom, G_CPLG_denom = sympy.fraction(G_CPLG)
        
        # export function for Denom and Nom
        G_CPLG_nom = sympy.lambdify(variable, G_CPLG_nom)
        G_CPLG_denom = sympy.lambdify(variable, G_CPLG_denom)
        
        # solving for stability
        solved_nom = sp.optimize.root(G_CPLG_nom, G_LS,tol=stab_tolerance)
        solved_denom = sp.optimize.root(G_CPLG_denom, G_LS, tol=stab_tolerance)

        # debug
        debug_info = "Delay: " + str(tau0) 
        
        # check zeros and plot
        nom_error = G_CPLG_nom(solved_nom.x)
        if solved_nom.success or (nom_error < stab_tolerance):
            zeros = solved_nom.x
        else:
            zeros = float("NaN")
            debug_info = debug_info + "\t Nom Error: " + str(nom_error)   
            debug_info = debug_info + "\t Sol: " + str(solved_nom.x)             
            #print(debug_info)
          
  
        # debug
        debug_info = "Delay: " + str(tau0) 
          
        # check poles and plot
        denom_error = G_CPLG_denom(solved_denom.x)
        if solved_denom.success or (denom_error < stab_tolerance):
            poles = solved_denom.x
        else:
            poles = float("NaN")
            debug_info = debug_info + "\t Denom Error: " + str(denom_error) 
            debug_info = debug_info + "\t Sol: " + str(solved_denom.x) 
            print(debug_info)
            
   
#############################################################################  
       
    # return variable
    next_state = []
    
    # Steady State Phase Model (2*pi*Freq) [x0, x1]
    next_state.append(omega0)
    next_state.append(omega1)

    # Steady State Phase Model (Delay) [tau0, tau1]
    next_state.append(tau0)
    next_state.append(tau1)
    
    if calc_stability:
        next_state.append(poles)
        next_state.append(zeros)
     
    # Triangular Function and their Derivertive
    next_state.append(delta_phase0)
    next_state.append(PD_Derivertive)
       
    return next_state

#############################################################################
###                 Fnet_Bi2N
#############################################################################
def Fnet_Bi2N_Stab(phase_delay, omega0_div, G_CPLG, variable, 
                   mode = 0, phaseshift=0, coupling_function= "triangle",
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
    phaseshift              (optionally) phase offset
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
    phase_error = mode - phase_delay - phaseshift*np.pi

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


#############################################################################
###                 Different Prototype Routines
#############################################################################
    
#############################################################################
###                 3th Generation, Homogen PCB V1.1
#############################################################################
    
def Model_3thGen_V1_1_homogen(phase_start=0, phase_end=2, phase_points=1e3,
                              VCO_freq0 = 2*np.pi*24.25e9,
                              Div_N=16*32,
                              K_VCO=2*np.pi*757.64e6,
                              K_PD = 1.6,
                              K_LF = 1,
                              coupling_function = "triangle",
                              phaseshift=-0.5,
                              invert=True,
                              laplace_var = ""):
    #############################################################################    
    """
   Simple Network of two coupled PLL System in Laplace Domain with default values from
   3. Gen Prototypes (Version 1.1) without the loopfilter by default.
    
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
    
    phaseshift              (optionally) additional phase shift
    invert                  (optionally) phase shift by pi
    laplace_var             (optionally) laplace variable (default: "s")
    """
    #############################################################################
    
    # laplace variable
    if not laplace_var:
        laplace_var = sympy.Symbol('s', positive=True)
    
      
    # feedforward and feedback transfer functions
    GTF_FF = K_PD * K_LF * (K_VCO/laplace_var)
    GTF_FB = 1/Div_N
     
    # open loop gain
    G_OL = K_PD * K_VCO *  1/Div_N
    print("Loop Gain: "  + str(G_OL/1e6) + "MHz and INV: " + str(invert))
    
    #############################################################################
    #           Two Coupled Nodes, Calculate Global Network Frequency
    #############################################################################
        
    # define phase range
    delay_phase = np.linspace(phase_start, phase_end, int(phase_points))
    
    # initial guess value
    freq_init = [VCO_freq0/Div_N, VCO_freq0/Div_N]
    
    #############################################################################
    
    # solving variable
    solved_freq = []
    solved_freq.append(freq_init)
    
    # solve iterative for each phase delay
    for delay in delay_phase:
        solved_freq.append(Fnet_Bi2N(solved_freq[-1],
                                     delay*np.pi,
                                     VCO_freq0/Div_N, 
                                     GTF_FF = GTF_FF,
                                     GTF_FB = GTF_FB,
                                     INV = invert, 
                                     phaseshift=phaseshift,
                                     variable = laplace_var, 
                                     coupling_function=coupling_function))
    
    # Output Format as Array/Matrix and delete init value
    solved_freq.pop(0)
    solved_freq = np.asarray(solved_freq)
    
    # extract data from node A and rescale to freq
    node_A_model = solved_freq[:,0]/(2*np.pi)
    node_A_tau_model = solved_freq[:,2]

    # poles and zeros
    poles = solved_freq[:,4]
    zeros = solved_freq[:,5]
 
    PD_CPLG = solved_freq[:,6]
    PD_CPLG_DT = solved_freq[:,7]

    #############################################################################
    

    #############################################################################    
    
    return [delay_phase, node_A_tau_model, node_A_model, poles, zeros, PD_CPLG, PD_CPLG_DT]
