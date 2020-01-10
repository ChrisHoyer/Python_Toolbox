#############################################################################
#   Different Scripts to Model State-Locked-Loops
#
#   - Solve_Models: Iterative solver for each SSL Model
#   - Bi_2Nodes_Lap_Freq: Simple Network in Laplace with two bidirectional coupled nodes to Calc Global Freq
#   - Bi_2Nodes_Tim_Freq: Script by Dimitrios, just for comparison
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Stand: 10.01.2020
#############################################################################

import numpy as np
import scipy as sp
#############################################################################
###                 Bi_2Nodes_Lap
#############################################################################
def Bi_2Nodes_Lap(last_state, delay, G_OL, delay_phase = True,
                  coupling_function= "sawthooth", 
                  coupling_scale=0.5, coupling_offset=1):
#############################################################################    
    """
   Simple Network in Laplace with two bidirectional coupled nodes to calculate
   the global network frequency
   Note, that this is a differential equation. Where the init value has to be set.
    
    parameters              description
    =====================  =============================================:
    last_state              last state from (as tuple!)
    delay                   current delay value
    G_OL                    open loop transfer function gain
    delay_phase             delay in phase in rad (N*pi)
    coupling_function       PD coupling function (XOR=sawthooth)
    coupling_scale          rescale of PD coupling function output
    coupling_offset         offset of PD coupling function output    
    
    return type
       next_state as tuple
    
    Example:
        
        import sll_toolbox as sll
        

    """
#############################################################################     
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
    phase_difference0 = last_state[0] + tau_phase0 - last_state[1]
    phase_difference1 = last_state[1] + tau_phase1 - last_state[0]
 

    if coupling_function == "sawthooth":
        
        # Calculate Sawthooth function        
        delta_phase0 = sp.signal.sawtooth(phase_difference0, width=0.5)
        delta_phase1 = sp.signal.sawtooth(phase_difference1, width=0.5)

        
    else:
        # no coupling
        delta_phase0 = phase_difference0
        delta_phase1 = phase_difference1
        
    # Scale coupling functions to your needs
    delta_phase0 = coupling_scale*(delta_phase0 + coupling_offset)
    delta_phase1 = coupling_scale*(delta_phase1 + coupling_offset)
    
    # return variable
    next_state = []
    
    # Steady State Frequency Model
    next_state[0] = G_OL * delta_phase0 + freq_0/N
    next_state[1] = G_OL * delta_phase1 + freq_0/N
        
    return next_state