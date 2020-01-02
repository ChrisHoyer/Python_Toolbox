#############################################################################
#
# - Extract_Sympy_1Var: Substitutes Sympy and generates numeric solution
# - BodePlot_FBCTRL: Generate BodePlot out of symbolic feedback transfer function
#
#   Autor: C. Hoyer (choyer.ch@gmail.com)
#   Stand: 02-01-2020
#############################################################################

import Basic_Toolbox as basic
import matplotlib.pyplot as plt
import numpy as np
import sympy

#############################################################################
#           Extract Magnitude and Phase from Laplace Transfer Function
#############################################################################
def Extract_Sympy_1Var(symbolic, eval_range, variable='s', evaluation="lambdify"):
############################################################################# 
    """
    Substitute Symbolic Function

    paramters              description
    =====================  =============================================:
    symbolic               Sympy function
    eval_range             subsitute content or evaluation range
    variable               (optional) variable parameter
    evaluation             (optional) define the used evalutaion function (lambdify | subs)
    
    return type
       array with values
       
    """   
 #############################################################################    
    # Generate Variable
    variable = sympy.Symbol(variable)
    
    # use Lambdify for fast results
    if evaluation == "lambdify":
        
        # Generate Function
        function = sympy.lambdify(variable, symbolic, modules="numpy")
    
        # Extract Function
        return function(eval_range)
    
    # use subs for precise results
    if evaluation == "subs":
        
        Solution = []
        
        # Evaluate:
        for eval_item in eval_range:
            evaluated = symbolic.evalf(subs={variable: eval_item})
            Solution.append(np.complex(evaluated))
       
        # return Values
        return Solution
            

#############################################################################
#           Generate BodePlot out of symbolic transfer function
#############################################################################
def BodePlot_FBCTRL(feedforward, feedback, freq, variable='s', evaluation="lambdify",
                    Add_LoopBW = True, Add_PhaseMargin = True,
                    Name_LoopBW= r'$f_\mathrm{loop}$ = ',Scale_LoopBW = 1e3,  Unit_LoopBW= 'kHz',
                    Name_PhaseMargin= r'$\varphi_\Delta$ = ', Max_PhaseMargin = -180,
                    Name_OL_dB = r'Open Loop Reference Phase $|\mathrm{G}_\mathrm{OL}|$',
                    Name_CL_dB = r'Closed Loop Reference Phase $|\mathrm{G}_\mathrm{CL}|$',
                    Name_OL_PH = r'Open Loop Reference Phase $\angle~(\mathrm{G}_\mathrm{OL})$',
                    Name_CL_PH = r'Closed Loop Reference Phase $\angle~(\mathrm{G}_\mathrm{CL})$'):
 
############################################################################# 
    """
    Generate BodePlot from Feedforward and Feedback Control
    
                           
    ------------->(+)-----> Feedforward  ------------->
                   ^ -                         |
                   |                           |
                   |                           |
                   -------- Feedback <----------
    

    paramters              description
    =====================  =============================================:
    feedforward             symbolic feedforward transfer function
    feedback                symbolic feedback transfer function
    freq                    frequency range of interest
    variable                (optional) laplace variable
    evaluation              (optional) define the used evalutaion function (lambdify|subs)  
    Add_LoopBW              (optional) Insert LoopBandwidth vertical line
    Name_LoopBW             (optional) Name of Loop Bandwidth
    Scale_LoopBW            (optional) Scale of Loop Bandwidth
    Unit_LoopBW             (optional) Unit of Loop Bandwidth
    Add_PhaseMargin         (optional) Insert Phase Margin horizontal line
    Name_PhaseMargin        (optional) Name of Phase Margin
    Max_PhaseMargin         (optional) Line at maximum Phase Margin
    Name_OL_dB              (optional) Name of Open Loop Magnitue (in dB)
    Name_CL_dB              (optional) Name of Closed Loop Magnitue (in dB)   
    Name_OL_PH              (optional) Name of Open Loop Phase (in degree)
    Name_CL_PH              (optional) Name of Closed Loop Phase (in degree)     
    
    return type
       array with values
       
    """   
 ############################################################################# 
    
    # Frequency Scale to j*w
    omega = 2j*np.pi * freq
    
    # Generate open and closed loop transfer function
    OpenLoop = feedforward * feedback
    ClosedLoop = feedforward / (1 + OpenLoop)
    
    # Generate OpenLoop Transfer Function
    OL_Extract = Extract_Sympy_1Var(OpenLoop, omega, evaluation=evaluation)
    OL_Extract = basic.CMPLX2Format(OL_Extract)
 
    # Generate ClosedLoop Transfer Function
    CL_Extract = Extract_Sympy_1Var(ClosedLoop, omega, evaluation=evaluation)
    CL_Extract = basic.CMPLX2Format(CL_Extract)
    
   
    # =================================== 
    # Plot Settings
    Xlabel_freq = ['', 'Hz']
    Ylabel_dB = ["", 'dB']
    Ylabel_phase = ["", '$^\circ$']   
    
    plot_mag = [[freq, OL_Extract['dB'], Name_OL_dB, 'linewidth=2.5'],
                [freq, CL_Extract['dB'], Name_CL_dB, 'linewidth=2.5']]
    
    plot_phase = [[freq, OL_Extract['PhaseDeg'], Name_OL_PH, 'linewidth=2.5'],
                  [freq, CL_Extract['PhaseDeg'], Name_CL_PH, 'linewidth=2.5']]

    # =================================== 
    # Generate Plot
    plt.figure(figsize=(10,10))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)   
    
    basic.SemiLogX_Plot(ax1, plot_mag, Xlabel_freq, Ylabel_dB)
    basic.SemiLogX_Plot(ax2, plot_phase, Xlabel_freq, Ylabel_phase)
    
    # =================================== 
    if Add_LoopBW | Add_PhaseMargin:
        
        # Find Gain of 1 in omega-domain and rescale
        loop_bw = sympy.solve(OpenLoop-1, variable)
        loop_bw = float(loop_bw[0])/(2*np.pi)
        
    # =================================== 
    if Add_LoopBW:
        
        # Rescale Loop Bandwidth
        loop_bw_str = loop_bw/Scale_LoopBW
        
        # Generate String
        loop_bw_str = str(np.round(loop_bw, decimals=3))
        
        # Add vertical line
        basic.Vline_Plot(ax1, loop_bw, Name_LoopBW + loop_bw_str + Unit_LoopBW)
        basic.Vline_Plot(ax2, loop_bw, Name_LoopBW + loop_bw_str + Unit_LoopBW)  
 
    # =================================== 
    if Add_PhaseMargin:

        # Find current phase
        current_phase = Extract_Sympy_1Var(ClosedLoop, (2j*np.pi*loop_bw))
        current_phase = np.angle(current_phase, deg=True)
        
        # calculate phase margin
        phase_margin = Max_PhaseMargin-current_phase
        
        # Add horizontal line
        basic.Hline_Plot(ax2, Max_PhaseMargin, str(Max_PhaseMargin) + '$^\circ$')
        basic.Hline_Plot(ax2, current_phase,  Name_PhaseMargin + str(phase_margin) + '$^\circ$', color='k', linestyle='--')

        
# =================================== 


    