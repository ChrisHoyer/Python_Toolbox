#############################################################################
#   Different Scripts to Model State-Locked-Loops
#
#   - Bi_2Nodes_Lap_Freq: Simple Network in Laplace with two bidirectional coupled nodes to Calc Global Freq
#   - Bi_2Nodes_Tim_Freq: Script by Dimitrios, just for comparison
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Stand: 10.01.2020
#############################################################################

import numpy as np
import sympy
import scipy as sp
from scipy import signal

import controltheory_toolbox as ctrl

#############################################################################
###                 Bi_2Nodes_Lap_Freq
#############################################################################
def Bi_2Nodes_Lap_Freq(last_state, delay, G_Gain, freq0_div,
                       INV = False, phaseshift=0, delay_phase = True, substitute = True,
                       calc_stability = True, GTF_FF = 0, GTF_FB = 0, variable='s',
                       coupling_function= "sawtooth", 
                       coupling_scale=0.5, coupling_offset=0):
#############################################################################    
    """
   Simple Network in Laplace with two bidirectional coupled nodes to calculate
   the global network frequency vs the phase delay
   Note, that this is a differential equation. Where the init value has to be set.
    
    parameters              description
    =====================  =============================================:
    last_state              last state from (as tuple!)
    delay                   current delay value
    GA_OL                   open loop transfer function (with laplace variable)
    freq0A_div              divided center frequency
    phaseshift              (optionally) phase offset related 
    delay_phase             (optionally) delay in phase in rad (N*pi)
    coupling_function       (optionally) PD coupling function (XOR=sawtooth, Mixer=cos)
    coupling_scale          (optionally) rescale of PD coupling function output
    coupling_offset         (optionally) offset of PD coupling function output 
    calc_stability          (optionally) calculate stability of the plots
    substitute              (optionally) sympy is better without numerics
    
    return type
       next_state and time delay as tuple [x0, x1, tau0, tau1]
    
    Example:
        
        import sll_toolbox as sll
        

    """
#############################################################################   

    # return variable
    next_state = []
    
    # delay format
    if delay_phase:
        # delay is respresented as phase in radian (0... 2pi)
        tau_phase0 = delay
        tau_phase1 = delay
    else:
        # convert delay into phase shift in radian at certain freq
        tau_phase0 = 2 * np.pi * last_state[0] * delay
        tau_phase1 = 2 * np.pi * last_state[1] * delay
    

    # Calculate phase difference based on delay in radian
    phase_difference0 = last_state[0] + tau_phase0 - last_state[1] - (INV*np.pi) - phaseshift
    phase_difference1 = last_state[1] + tau_phase1 - last_state[0] - (INV*np.pi) - phaseshift
 

    if coupling_function == "sawtooth":
        
        # Calculate Sawthooth function        
        delta_phase0 = sp.signal.sawtooth(phase_difference0, width=0.5)
        delta_phase1 = sp.signal.sawtooth(phase_difference1, width=0.5)

    elif coupling_function == "cos":
        
        # Calculate Sawthooth function        
        delta_phase0 = np.cos(phase_difference0)
        delta_phase1 = np.cos(phase_difference1)
        
    else:
        # no coupling
        delta_phase0 = phase_difference0
        delta_phase1 = phase_difference1
        
    # Scale coupling functions to your needs
    delta_phase0 = coupling_scale*(delta_phase0 + coupling_offset)
    delta_phase1 = coupling_scale*(delta_phase1 + coupling_offset)
    
    # calculate global frequency
    freq0 = G_Gain * delta_phase0 + freq0_div
    freq1 = G_Gain * delta_phase1 + freq0_div
  
    # calculate stability of network frequency
    if calc_stability:
        
        # phase difference, sign
        ddt = np.sign(freq0 - last_state[0])
        
        # Generate transfer function dependcy on phase difference
        GTF_FF_Delay = GTF_FF * ddt
        G_TF = GTF_FF_Delay / (1+ GTF_FF_Delay * GTF_FB)
        
        # substitue
        [G_TF, Mapping] = ctrl.Substitute_Datatype(G_TF, datatype="Float")   
    
        # Calculate Poles and Zeros
        G_TF = sympy.simplify(G_TF)
        nom, denom = sympy.fraction(G_TF)
        poles = sympy.roots(denom, variable) 
        zeros = sympy.roots(nom, variable)
        
        # Resubstitute and rearrange poles
        value = []
        for key in poles.keys():
            resub = ctrl.ReSubstitute_Datatype(key, Mapping)
            # convert to complex number
            value.append(complex(resub.evalf()))  
        poles = value
 
        # Resubstitute and rearrange poles
        value = []
        for key in zeros.keys():
            resub = ctrl.ReSubstitute_Datatype(key, Mapping)
            # convert to complex number
            value.append(complex(resub.evalf()))  
        zeros = value       
        
    # calculate time delay
    tau0 = delay/(2*np.pi*freq0)
    tau1 = delay/(2*np.pi*freq1)
    
    # Steady State Phase Model (Freq) [x0, x1]
    next_state.append(freq0)
    next_state.append(freq1)

    # Steady State Phase Model (Delay) [tau0, tau1]
    next_state.append(tau0)
    next_state.append(tau1)
    
    if calc_stability:
        next_state.append(poles)
        next_state.append(zeros)
            
    return next_state

#############################################################################  