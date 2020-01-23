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
import scipy as sp
from scipy import signal

#############################################################################
###                 Bi_2Nodes_Lap_Freq
#############################################################################
def Bi_2Nodes_Lap_Freq(last_state, delay, GA_OL, freq0A_div,
                       GB_OL = 0.0, freq0B_div = 0.0, INV = True, shift=0,
                       delay_phase = True, coupling_function= "sawthooth", 
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
    GA_OL                   open loop transfer function gain node A
    GB_OL                   (optionally) open loop transfer function gain node B
    freq0A_div              divided center frequency node A
    freq0B_div              (optionally) divided center frequency node B
    phaseoffset_A           (optionally) phase offset related to node A
    phaseoffset_B           (optionally) phase offset related to node B
    delay_phase             (optionally) delay in phase in rad (N*pi)
    coupling_function       (optionally) PD coupling function (XOR=sawthooth, Mixer=cos)
    coupling_scale          (optionally) rescale of PD coupling function output
    coupling_offset         (optionally) offset of PD coupling function output    
    
    return type
       next_state and time delay as tuple [x0, x1, tau0, tau1]
    
    Example:
        
        import sll_toolbox as sll
        

    """
#############################################################################   
    
    # heterogeneous or homogeneous case
    if GB_OL == 0.0:
        GB_OL = GA_OL
    
    if freq0B_div == 0.0:
        freq0B_div = freq0A_div
        
    # Consider a frequency difference
    delta_freq0 = np.abs(freq0B_div - freq0A_div)
    mean_freq0 = np.mean([freq0B_div,freq0A_div])
    phaselag = 0.25*np.pi
    
    
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
    phase_difference0 = last_state[0] - tau_phase0 - last_state[1] - (INV*np.pi) - shift
    phase_difference1 = last_state[1] - tau_phase1 - last_state[0] - (INV*np.pi) - shift
 

    if coupling_function == "sawthooth":
        
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
    freq0 = GA_OL * delta_phase0 + freq0A_div
    freq1 = GB_OL * delta_phase1 + freq0B_div
    
    # calculate time delay
    tau0 = delay/(2*np.pi*freq0)
    tau1 = delay/(2*np.pi*freq1) 
    
    # return variable
    next_state = []
    
    # Steady State Phase Model (Freq) [x0, x1]
    next_state.append(freq0)
    next_state.append(freq1)

    # Steady State Phase Model (Delay) [tau0, tau1]
    next_state.append(tau0)
    next_state.append(tau1)
            
    return next_state

#############################################################################  