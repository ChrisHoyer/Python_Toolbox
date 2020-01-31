#############################################################################
#   Different Scripts to Model State-Locked-Loops
#
#   - Bi_2Nodes_Freq_Stab: Simple Network in Laplace with two bidirectional coupled nodes to Calc Global Freq
#
#   Prototype Models
#   - Model_3thGen_V1_1_homogen
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Stand: 10.01.2020
#############################################################################
import controltheory_toolbox as ctrl
import basic_toolbox as basic

import numpy as np
import sympy

import scipy as sp
from scipy import signal

import matplotlib.pyplot as plt


#############################################################################
###                 Bi_2Nodes_Freq_Stab
#############################################################################
def Bi_2Nodes_Freq_Stab(last_state, delay, G_Gain, freq0_div,
                       INV = False, phaseshift=0, delay_phase = True, substitute = True,
                       calc_stability = True, GTF_FF = 0, GTF_FB = 0, variable='s',
                       Invert_FF = False, coupling_function= "sawtooth", 
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
    phase_difference0 = last_state[0] + tau_phase0 - last_state[1]
    phase_difference1 = last_state[1] + tau_phase1 - last_state[0]
 
    # Add Invertion and Phaseshifts (Represents varphi_const), which is equal
    phase_difference0 = phase_difference0  - (INV*np.pi) - phaseshift
    phase_difference1 = phase_difference1  - (INV*np.pi) - phaseshift

    # coupling function, XOR
    if coupling_function == "sawtooth":
        
        # Calculate Sawthooth function        
        delta_phase0 = sp.signal.sawtooth(phase_difference0, width=0.5)
        delta_phase1 = sp.signal.sawtooth(phase_difference1, width=0.5)

    # coupling function, MIXER
    elif coupling_function == "cos":
        
        # Calculate Sawthooth function        
        delta_phase0 = np.cos(phase_difference0)
        delta_phase1 = np.cos(phase_difference1)
    
    # no coupling function, XOR    
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
        
        # phase shift or inversion in feed forward
        if Invert_FF:
            ddt = ddt * -1
           
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
    

#############################################################################
###                 Different Prototype Routines
#############################################################################
    
#############################################################################
###                 3th Generation, Homogen PCB V1.1
#############################################################################
    
def Model_3thGen_V1_1_homogen(phase_start=0, phase_end=2, phase_points=1e3,
                              Div_N=16*32, VCO_freq0 = 24.25e9,  invert=False,
                              invert_ff=True, plot=False):
    #############################################################################
    #           System Setup - Model
    #############################################################################
    # laplace variable
    s = sympy.Symbol('s', positive=True)
    
    # Sensitivity of the VCO [Hz/V]
    K_VCO = 757e6
    K_VCO = 3e9*1/(2*np.pi)
    
    # Sensitivity of the PD [V/rad]
    K_PD = 0.8/(1*np.pi)
    K_PD = 1
    
    # Divider
    N = Div_N
    
    # Loop Filter
    f_LF = 1e6
    K_LF = 1/( 1 + (s * 1/(2*np.pi*f_LF) ))
    K_LF = 1
        
    # open loop gain
    G_OL = K_PD * 2 * np.pi * K_VCO *  1/N
    
    # feedforward and feedback transfer functions
    GTF_FF = K_PD * K_LF * 2 * np.pi * K_VCO * 1/s
    GTF_FB = 1/N
    
    # Invert FF (DiffSe: 180, Adder: 180, Offset:180)
    #invert_ff = True

    # Center frequency of the System
    freq_0 = VCO_freq0
    
    #print loop gain
    print("Loop Gain: "  + str(G_OL/1e6) + "MHz/pi")
    
    #############################################################################
    #           Two Coupled Nodes, Calculate Global Network Frequency
    #############################################################################
        
    # define phase range
    delay_phase = np.linspace(phase_start, phase_end, phase_points)
    
    # initial value
    freq_init = [freq_0/N, freq_0/N]
    
    #############################################################################
    
    # solving variable
    freq_out = []
    freq_out.append(freq_init)
    
    # solve iterative
    for delay in delay_phase:
        freq_out.append(Bi_2Nodes_Freq_Stab(freq_out[-1], delay*np.pi,
                                           G_OL, freq_0/N, 
                                           GTF_FF = GTF_FF, GTF_FB=GTF_FB,
                                           INV = invert,
                                           variable=s, Invert_FF=invert_ff, 
                                           coupling_function="sawtooth"))
    
    # Output Format as Array/Matrix and delete init value
    freq_out.pop(0)
    freq_out = np.asarray(freq_out)
    
    # extract data from node A
    node_A_model = freq_out[:,0]
    node_A_tau_model = freq_out[:,2]

    # find max. pol Value
    max_real_pol = [np.max(np.real(item)) for item in list(freq_out[:,4])]

    #############################################################################
    #           Plot Global Frequency
    #############################################################################
    if plot:    
    
        # plot results
        Xlabel_time = [r'Delay $\tau_\mathrm{AB}$', 's']
        Xlabel_phase = [r'Phase Difference $\varphi_\mathrm{AB}$', r'$\pi$']
    
        Ylabel = [r'Network Frequency $f_\mathrm{NET}$', 'Hz']
    
        plot1 = [[delay_phase, node_A_model, r'Node A ($f_\mathrm{OUT}^\mathrm{A}$)']]
        
        plot1A = [[delay_phase, max_real_pol, r'Node A Pol', 'linestyle=--']]
    
        plot2 = [[node_A_tau_model, node_A_model, r'Node A ($f_\mathrm{OUT}^\mathrm{A}$)']]
    
    
        # Generate Plot
        plt.figure(figsize=(10,8))
        ax1 = plt.subplot(311)
        ax1a = ax1.twinx() 
        
        basic.Linear_Plot(ax1, plot1, Xlabel_phase, Ylabel, LegendLoc=0) 
        basic.Linear_Plot(ax1a, plot1A, Xlabel_phase, Ylabel,TwinX=ax1, LegendLoc=0)  
        
        ax2 = plt.subplot(312)
        basic.Linear_Plot(ax2, plot2, Xlabel_time, Ylabel, LegendLoc=0)

    #############################################################################    
    
    return [delay_phase, node_A_tau_model, node_A_model, max_real_pol]