#############################################################################
#
# - Extract_Sympy_1Var: Substitutes Sympy and generates numeric solution
# - BodePlot_FBCTRL: Generate BodePlot out of symbolic feedback transfer function
# - StepResponse: Generate Step Response with Heaviside Fct from symbolic transfer function
# - Substitute_Datatype: Substitute constant values with symboles
# - ReSubstitute_Datatype: Resubstitute constant values with symboles
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Stand: 03-01-2020
#############################################################################

import basic_toolbox as basic
import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy.integrals import inverse_laplace_transform
from sympy.parsing.sympy_parser import parse_expr

#############################################################################
#           Extract Magnitude and Phase from Laplace Transfer Function
#############################################################################
def Extract_Sympy_1Var(symbolic, eval_range, variable='s', evaluation="lambdify",
                       modules=['numpy', 'sympy'], Heaviside_Value=0):
############################################################################# 
    """
    Substitute Symbolic Function

    paramters              description
    =====================  =============================================:
    symbolic               Sympy function
    eval_range             subsitute content or evaluation range
    variable               (optional) variable parameter
    evaluation             (optional) define the used evalutaion function (lambdify | subs)
    Heaviside_Value        (optional) Value for Heaviside Function at H(0)
    
    return type
       array with values
       
    """   
 #############################################################################   
    def SetHeaviside(symbolic, Heaviside_Value=0):
        
        # Element is a Heaviside Element?
        if str(type(symbolic)) == "Heaviside":
            # Extract old Argument and generate new one
            symbolic_args = (symbolic.args[0], Heaviside_Value)
            symbolic = symbolic.func(*symbolic_args)
            # return back into parent
            return symbolic
            
        # Save all return args and write back into parent
        new_args = ()
    
        # Iterate the tree
        for arg in symbolic.args:
                 
            # Argument has another arguments?
            if arg.args:
                arg = SetHeaviside(arg, Heaviside_Value=Heaviside_Value)
            
            # write back into new arguments
            new_args = new_args + (arg,)

        # Write new arguments back into parent
        symbolic = symbolic.func(*new_args)
                                    
        return symbolic

    
 #############################################################################  
    # Generate Variable
    if isinstance(variable, str):
        variable = sympy.Symbol(variable)
        
        
    # Set Heaviside to default value at jump
    symbolic = SetHeaviside(symbolic)
                    
    # use Lambdify for fast results
    if evaluation == "lambdify":
        
        # Generate Function
        function = sympy.lambdify(variable, symbolic, modules=modules)
    
        try:
            # use array directly
            output = function(eval_range)
            
        except:
            # iterate each item
            output = []
            for eval_item in eval_range:
                evaluated = function(eval_item)  
                output.append(float(evaluated))
                
        # Extract Function
        return output
    
    # use subs for precise results
    if evaluation == "subs":
        
        Solution = []
        
        # Evaluate:
        for eval_item in eval_range:
            evaluated = symbolic.evalf(subs={variable: eval_item})
            evaluated = np.complex(evaluated)
       
            Solution.append(evaluated)
       
        # return Values
        return Solution
            
#############################################################################
#           Generate BodePlot out of symbolic transfer function
#############################################################################
def BodePlot_FBCTRL(feedforward, feedback, freq, variable='s', evaluation="lambdify",
                    Add_LoopBW = False, Add_PhaseMargin = True, Add_MaxPhaseMargin = True,
                    PlotBlackWhite = False,
                    Name_LoopBW= r'$f_\mathrm{loop}$ = ',Scale_LoopBW = 1e3,  Unit_LoopBW= 'kHz',
                    Name_PhaseMargin= r'$\varphi_\Delta$ = ', Max_PhaseMargin = -180,
                    Name_OL_dB = r'Open Loop Reference Phase $|\mathrm{G}_\mathrm{OL}|$',
                    Name_CL_dB = r'Closed Loop Reference Phase $|\mathrm{G}_\mathrm{CL}|$',
                    Name_OL_PH = r'Open Loop Reference Phase $\angle~(\mathrm{G}_\mathrm{OL})$',
                    Name_CL_PH = r'Closed Loop Reference Phase $\angle~(\mathrm{G}_\mathrm{CL})$',
                    Xlabel_freq = ['', 'Hz'], Ylabel_dB = ["", 'dB'], Ylabel_PH = ["", '$^\circ$'] ):
 
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
    PlotBlackWhite          (optional) plot only in black and white
    evaluation              (optional) define the used evalutaion function (lambdify|subs)  
    Add_LoopBW              (optional) Insert LoopBandwidth vertical line
    Name_LoopBW             (optional) Name of Loop Bandwidth
    Scale_LoopBW            (optional) Scale of Loop Bandwidth
    Unit_LoopBW             (optional) Unit of Loop Bandwidth
    Add_PhaseMargin         (optional) Insert Phase Margin horizontal line
    Add_MaxPhaseMargin      (optional) Insert Phase Margin horizontal line at Max_PhaseMargin
    Name_PhaseMargin        (optional) Name of Phase Margin
    Max_PhaseMargin         (optional) Line at maximum Phase Margin
    Name_OL_dB              (optional) Name of Open Loop Magnitue (in dB)
    Name_CL_dB              (optional) Name of Closed Loop Magnitue (in dB)   
    Name_OL_PH              (optional) Name of Open Loop Phase (in degree)
    Name_CL_PH              (optional) Name of Closed Loop Phase (in degree)     
    Xlabel_freq             (optional) Label and Unit of Frequency X-Axis
    Ylabel_dB               (optional) Label and Unit of dB Y-Axis 
    Ylabel_PH               (optional) Label and Unit of phase Y-Axis  
    
    return type
       plot
       
    """   
 ############################################################################# 
    # Generate Variable
    variable = sympy.Symbol(variable)
    
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
    
    plot_mag = [[freq, OL_Extract['dB'], Name_OL_dB, 'linewidth=2.5'],
                [freq, CL_Extract['dB'], Name_CL_dB, 'linewidth=2.5']]
    
    plot_phase = [[freq, OL_Extract['PhaseDeg'], Name_OL_PH, 'linewidth=2.5'],
                  [freq, CL_Extract['PhaseDeg'], Name_CL_PH, 'linewidth=2.5']]

    # =================================== 
    # Generate Plot
    plt.figure(figsize=(10,10))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)   
    
    basic.SemiLogX_Plot(ax1, plot_mag, Xlabel_freq, Ylabel_dB, BlackWhite=PlotBlackWhite)
    basic.SemiLogX_Plot(ax2, plot_phase, Xlabel_freq, Ylabel_PH, BlackWhite=PlotBlackWhite)
    
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
        loop_bw_str = str(np.round(loop_bw_str, decimals=3))
        
        if PlotBlackWhite:
            # Add vertical line
            basic.Vline_Plot(ax1, loop_bw, Name_LoopBW + loop_bw_str + Unit_LoopBW, color="gray", linestyle=(0, (5, 10)))
            basic.Vline_Plot(ax2, loop_bw, Name_LoopBW + loop_bw_str + Unit_LoopBW, color="gray", linestyle=(0, (5, 10)))  
        else:
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
        phase_margin_str = str(np.round(phase_margin, decimals=3))
        
        if PlotBlackWhite:
            # Add horizontal line
            if Add_MaxPhaseMargin:
                basic.Hline_Plot(ax2, Max_PhaseMargin, str(Max_PhaseMargin) + '$^\circ$', color="gray", linestyle=(0, (5, 10)))
            basic.Hline_Plot(ax2, current_phase,  Name_PhaseMargin + phase_margin_str + '$^\circ$', color='gray', linestyle=(0, (5, 10)))

        else:
            # Add horizontal line
            if Add_MaxPhaseMargin:
                basic.Hline_Plot(ax2, Max_PhaseMargin, str(Max_PhaseMargin) + '$^\circ$')
            basic.Hline_Plot(ax2, current_phase,  Name_PhaseMargin + phase_margin_str + '$^\circ$', color='k', linestyle='--')
            

   
    # return plot
    return plt

#############################################################################
#           Generate BodePlot out of symbolic transfer function
#############################################################################
def BodePlot(systems, evaluation="lambdify", PlotBlackWhite=False, plotphase=True,
             Xlabel_freq = ['', 'Hz'], Ylabel_dB = ["", 'dB'], Ylabel_PH = ["", '$^\circ$'],
             fig = plt.figure(figsize=(10,10))):
 
############################################################################# 
    """
    Generate BodePlot from transfer function
    

    paramters              description
    =====================  =============================================:
    systems                 symbolic systems transfer function as List
    evaluation              (optional) define the used evalutaion function (lambdify|subs) 
    PlotBlackWhite          (optional) plot only in black and white
    Xlabel_freq             (optional) Label and Unit of Frequency X-Axis
    Ylabel_dB               (optional) Label and Unit of dB Y-Axis 
    Ylabel_PH               (optional) Label and Unit of phase Y-Axis
    fig                     (optional) matplot figure
    
    return type
       plot
       
    """   
 ############################################################################# 

    # Plot Settings
    plot_mag = []
    plot_phase = []

    # ===================================     
    # Iterate all Equations
    for index in range(len(systems)):
        
        # get system
        system = systems[index]
        
        # Frequency Scale to j*w
        freq = system[0]
        omega = 2j*np.pi * freq
        
        # Generate OpenLoop Transfer Function
        system_extract = Extract_Sympy_1Var(system[1], omega, evaluation=evaluation)
        system_extract = basic.CMPLX2Format(system_extract)
 
        # Generate Magnitude Plot Settings
        plot_mag.append([freq, 
                         system_extract['dB'],
                         system[2] + " " + r'$|' + system[3] + r'|$',
                         system[4]
                         ])
        
        if plotphase:
            # Generate Magnitude Plot Settings
            plot_phase.append([freq, 
                               system_extract['PhaseDeg'],
                               system[2] + " " + r'$\angle~(' + system[3] + r'(|)$',
                               system[4]
                               ]) 
   
    # =================================== 
    # Generate Plot
    
    # plot with phase
    if plotphase:
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)
        fig.add_axes(ax1)
        fig.add_axes(ax2)
         
        basic.SemiLogX_Plot(ax1, plot_mag, Xlabel_freq, Ylabel_dB, BlackWhite=PlotBlackWhite)
        basic.SemiLogX_Plot(ax2, plot_phase, Xlabel_freq, Ylabel_PH, BlackWhite=PlotBlackWhite)
            
    else:
        ax1 = plt.subplot(111)
        fig.add_axes(ax1)
        
        
        basic.SemiLogX_Plot(ax1, plot_mag, Xlabel_freq, Ylabel_dB, BlackWhite=PlotBlackWhite)      
       
    # return plot
    return plt
 
     
#############################################################################
#           Generate StepResponse out of symbolic transfer function
#############################################################################
def StepResponse(system, time, delay=1, variable_laplace='s', variable_time='t',
                 substitute=True, plot=True, normalize_plot=True,
                 evaluation="subs", Heaviside_Value=0,
                 Name_StepFct = r'Step Function, normalized',
                 Name_System = r'Step Response, normalized',
                 Xlabel_time = ['', 's'], Ylabel_Amplitude = ["", 'V']  
                 ):
    
############################################################################# 
    """
    Generate Step Response from System
    

    paramters              description
    =====================  =============================================:
    system                  system symbolic transfer function
    time                    time array for evaluation
    delay                   (optional) delay of heaviside function
    variable_laplace        (optional) variable in laplace domain
    variable_time           (optional) variable in time domain
    substitute              (optional) substitute floats during inverse laplace
    plot                    (optional) plot step response output
    normalize_plot          (optional) normalize end value to maximum
    evaluation              (optional) define the used evalutaion function (lambdify|subs)  
    Heaviside_Value         (optional) Value for Heaviside Function at H(0)
    Name_StepFct            (optional) Name of StepFunction Label
    Name_System             (optional) Name of System response label    
    Xlabel_time             (optional) Label and Unit of Time X-Axis
    Ylabel_Amplitude        (optional) Label and Unit of Amplitude Y-Axis 
    
    
    return type
       symbolic time domain functions
       
    """   
 ############################################################################# 
    
    # ===================================     
    # Generate symbolic variables (time only positive)
    variable_laplace = sympy.Symbol(variable_laplace)
    variable_time = sympy.Symbol(variable_time, positive=True)
    
    # Generate Heaviside function with delay
    heaviside_laplace = 1/variable_laplace * sympy.exp(-delay*variable_laplace)
    heaviside_time = inverse_laplace_transform(heaviside_laplace, variable_laplace, variable_time).doit()
    
    # =================================== 
    # Step Response in Laplace Domain
    StepReponse_laplace = system * heaviside_laplace
    
    # Substitute all float numbers
    if substitute:
         [StepReponse_laplace, Mapping] = Substitute_Datatype(StepReponse_laplace, datatype="Float")   
    
    # Convert StepResponse into Time Domain
    StepReponse_time = inverse_laplace_transform(StepReponse_laplace, variable_laplace, variable_time).doit()
    
    # Re-Substitute all float numbers
    if substitute:
         StepReponse_time = ReSubstitute_Datatype(StepReponse_time, Mapping)   
    
    # =================================== 
    if plot:   

        # evaluation of that data
        heaviside = Extract_Sympy_1Var(heaviside_time, time, variable=variable_time,
                                       evaluation=evaluation,  Heaviside_Value=Heaviside_Value)   
        stepresponse = Extract_Sympy_1Var(StepReponse_time, time, variable=variable_time,
                                          evaluation=evaluation, Heaviside_Value=Heaviside_Value) 
        
        # As Array
        heaviside = np.asarray(heaviside)
        stepresponse = np.asarray(stepresponse)
         
        # normalize
        if normalize_plot:
            heaviside = heaviside/heaviside[-1]
            stepresponse = stepresponse/stepresponse[-1]
              
    # ===================================        
        # Plot Settings
    
        plot_amp = [[time, heaviside, Name_StepFct, 'linewidth=0.5, marker=x, markersize=10'],
                   [time, stepresponse, Name_System, 'linewidth=1']]

    # =================================== 
        # Generate Plot
        plt.figure(figsize=(10,10))
        ax1 = plt.subplot(111)
    
        basic.Linear_Plot(ax1, plot_amp, Xlabel_time, Ylabel_Amplitude)
        
        plt.show()
    
        
    # =================================== 
    # return time domain systems
    return[heaviside_time, StepReponse_time]

#############################################################################
#           Substitute constant values with symboles
#############################################################################     
def Substitute_Datatype(symbolic_expression, datatype="Float", add_sym="real=True"):
        
    # export expression as text
    expression = sympy.srepr(symbolic_expression)
        
    # Extract Datatype-Command
    splitted_Expression = expression.split(datatype)
    splitted_Expression.pop(0)
        
    # All found commands
    mapping = {}
    
    # Generate Output Type
    expression_new = expression
    
    # Datatype detected?
    if len(splitted_Expression) > 0:
        
        # Generate Mapping
        for i in range(len(splitted_Expression)):
            
            # Extract Datatype
            content = splitted_Expression[i][splitted_Expression[i].find("(")+
                                                 1:splitted_Expression[i].find(")")]
            content = (datatype + "(" + content + ")")
            
            # Generate new Symbol string
            symbol = "Symbol('" + str(datatype[0]) + str(i) + "', "+ add_sym + ")"
                
            # add to mapping
            mapping[symbol] = content
                
            # Replace in string
            expression_new = expression_new.replace(content, symbol)
        
        # import expression as nummeric equation
        parsed = parse_expr(expression_new) 
            
    else:
            
        # return symbolic_expression nothing found
        parsed = symbolic_expression
         
    # return new expression and mapping
    return [parsed, mapping]
    
#############################################################################
#           Resubstitute constant values with symboles
#############################################################################   
def ReSubstitute_Datatype(symbolic_expression, mapping):
        
    # export expression as text
    expression = sympy.srepr(symbolic_expression)
 
    if len(mapping) > 0:
        
        # iterate all keys
        for mapping_key in mapping.keys():
            
            # Replace in string
            expression = expression.replace(mapping_key, mapping[mapping_key])
            
        # import expression as nummeric equation
        parsed = parse_expr(expression) 
        
    else:
        # return symbolic_expression nothing found
        parsed = symbolic_expression
    
    # return old expression
    return parsed           
