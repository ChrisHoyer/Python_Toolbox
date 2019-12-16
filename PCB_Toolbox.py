#############################################################################
#   1.  CSV2Dict - script for importing CSV Files
#   2.  CSV2Area - script for calculating polygon area based on coordinates
#   3.  FitFct_Exp - script that fits data against an Exp Function
#   4.  FindPoint_FitFct - fits higher order Polynom against Dataset to find point
#   5.  Linearization_Point - linearization around one point
#   6.  Linear_Plot - linear plot function with automated labeling
#   7.  SemiLogX_Plot - semilog x plot function with automated labeling
#	8.  Vline_Plot - generates vertical line in plot with label
#	9.  Hline_Plot - generates horizontal line in plot with label
#   10. Align_YAxis - Align two YGrids in one plot
#	11. FindPoint_NextValue - find nearest point
#   12. Digitalize_Data - Generates binary stream from data and clock
#	13. CMPLX2Format - converts complex numbers in different formats
#   14. Average - Average a Dataset
#   15. FFTData - Generates a FFT out the Data
#   16. String2List - Generates a List out of a csv with one line
#
#   Autor: C. Hoyer (choyer.ch@gmail.com)
#   Stand: 09-12-2019
#############################################################################

from scipy.optimize import curve_fit
import matplotlib.ticker as tck
import matplotlib.pyplot as plt
import numpy as np
import csv
import re

#############################################################################
###         Import CSV File to a dictionary
#############################################################################
def CSV2Dict(csv_file, delimiter=';', complexdelimiter='/', 
             headerline_ends= 1, blockheader=False, dictkey=False,
             headerkeys=[], **kwargs):
############################################################################# 
        """
    Imports all CSV Data (also complex and logarithmic data)
    Specialized for LTspice and ADS

    paramters              description
    =====================  =============================================:
    csv_file                csv file
    delimiter               (option) csv file delimiter
    complexdelimiter        (option) delimiter if there is a complex value
    blockheader             (option) each block has a header? (e.g. ADS, or global header?)
    headerline_ends         (option) end of header (skip until this line)
    dictkey                 (option) try to parse dictkeys instead of block number
    headerkeys              (option) custom header names (must be equal to row count)
    
    return type
       dictionary with each column or nested dictionary
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # import csv
        file = 'test.csv'
        
        # calculate
        Data = pcb.CSV2Dict(file)
         
   
    """
#############################################################################  
        # row state
        States = ['', 'SimpleNumber', 'ComplexNumber', 'NewBlock']
        
        # Matching Numbers to parse float / scientific data
        scientific_number = re.compile('[-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?')  

############################################################################# 
        
        def CSVData_DefineContent(header, row, v):
            
            # Try to Parse the Data
            Current_State = States[0]
            
            # only plain Number (normal or scientific)
            try:
                float(v)
                Current_State = States[1]
                    
            except:
                
                # Row without Data?
                if len(row) != len(header):
                    Current_State = States[3]  
                    
                # Row with Text?
                if (type(v) is str) & (v != "") :
                     Current_State = States[3]                     
                    
                # Contains two values? -> Complex
                elif len(v.split(complexdelimiter)) == 2:
                    Current_State = States[2]

            return Current_State

############################################################################# 
            
        def CSVData_ParseContent(header, row, v):
            
            # Classify Data
            State = CSVData_DefineContent(header, row, v)
                        
            try:
                # Convert Float Data
                if State == States[1]:
                    # Only a number?:
                    v_parsed = float(v)
                
                # Convert Complex Data
                elif State == States[2]:
                    # Split
                    v = v.split(complexdelimiter)
                        
                    # check length 
                    if len(v) == 2:
                            
                        # convert to magniture and angle
                        v_mag = [float(x) for x in re.findall(scientific_number, v[0])][0]
                        v_angl =[float(x) for x in re.findall(scientific_number, v[1])][0]
                            
                        # check if in dB or linear and convert to volts
                        if v[0].find("dB", re.IGNORECASE):
                            v_mag = 10**(v_mag/20)
                                
                        # polar coordinates to cartesian
                        v_parsed = v_mag * np.exp(1j*v_angl)
                                
                # Convert Header Data                                
                elif State == States[3]:
                    # new Block starts
                    v_parsed = States[3]
                    
                else:
                    v_parsed = float('NaN')
                    print("Data " + str(v) + " could not be parsed!")                  
                
            except:
                print("Data " + str(v) + " could not be parsed!") 
                
            return v_parsed
            
#############################################################################   
        # Import CSV Matrix
        csv_reader = csv.reader(open(csv_file), delimiter=delimiter)
        global_line_count = 0
        
        # generate Block
        Block = {}
        line_count = 0
        NewBlockStart = False
        
        # generate Block Content
        block_name = ''
        block_count = 0
        Data = {}
        
        # no data region
        data_region = False
        
 #############################################################################  
        # Import Block
        for row in csv_reader:
              
            # empty row? -> skip and  new Block
            if not row:
                NewBlockStart = True
                continue
            
            # get Block
            if NewBlockStart:
                
                # Blockname found?
                if not dictkey:
                    block_name = block_count
                   
                # delete rows
                block_deleteRow = []
                    
                # check if rows are empty:
                for block_row in Block:
                    if len(Block[block_row]) == 0:
                        block_deleteRow.append(block_row)
                    
                # delete empty rows
                for deleteRow in block_deleteRow:
                    Block.pop(deleteRow)
                        
                # add if it is not empty
                if Block:              
                    # Add to Data
                    Data[block_name] = Block
                
                # clear flag
                NewBlockStart = False

                # start new lines and count Block
                Block = {}
                block_name = ''
                block_count += 1
                line_count = 0
            

            # Skip Global Header
            if  global_line_count < (headerline_ends):
                global_line_count += 1
                continue

            # Find beginning of the Block
            if line_count == 0 :
                
                if blockheader or block_count == 0:
                    # read new/first block header
                    header = [r.replace(' ', '') for r in row]
                
                # Generate Custom Header
                if len(header) == len(headerkeys):
                      header = headerkeys
                
                # Generate Columnheader from first row
                for h in header:
                    Block[h] = [] 
            
            # decide when the data region beginns
            if not(blockheader):
                data_region = True
            else:
                data_region = False

            # end of header reached? -> Data region
            if data_region:
                
                # iterate row
                for h, v in zip(header, row):
                    
                    # Classify Data and Parse
                    parsed_Data = CSVData_ParseContent(header, row, v)
                    
                    # Header found?
                    if parsed_Data == States[3]:
                        NewBlockStart = True
                        continue
                        
                    Block[h].append(parsed_Data)
                              
            # next line    
            line_count += 1
            global_line_count += 1
            
        # last Block
        if Block:
            
            # Blockname?
            if not dictkey:
                block_name = block_count
            
            # Add to Dataset
            Data[block_name] = Block
          
        # only one Block? 
        if block_count == 0:
            return Block
        else:
            return Data
############################################################################# 
    
#############################################################################
###         Import CSV File and Calculate
#############################################################################
def CSV2Area(coordinate_file, draw_polygon=False, delimiter=';', **kwargs):
#############################################################################    
    """
    Calculates an Area given by CSV coordinates

    paramters              description
    =====================  =============================================:
    coordinate_file         csv file
    draw_polygon            (option) plot the input polygon
    delimiter               (option) csv file delimiter
    
    return type
       float with the calculated area
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # import csv
        file = 'test.csv'
        
        # calculate and print
        Area = pcb.CSV2Area(file, draw_polygon=True)
        
        print('This Polygon has an Area of ' + str(Area) + ' mm2')      
   
    """
#############################################################################         

    # import Data
    csv_data = CSV2Dict(coordinate_file)
       
    # X and Y Coordinate List     
    X = csv_data['X']
    Y = csv_data['Y']
            
    # draw Polygon
    if draw_polygon:
        
        # generate closed shape
        x_draw = X.copy()
        x_draw.append(X[0])
        y_draw = Y.copy()
        y_draw.append(Y[0])  
        
        # simple plot
        plt.figure()
        plt.plot(x_draw, y_draw)
        plt.show
        
    # calculate area with Shoelace formula (Trapz in Matlab)
    
    # move index of Y+1, X+1 respectively
    Roll_Y = np.dot( X, np.roll(Y,1))
    Roll_X = np.dot(Y, np.roll(X,1))
    Area = np.abs(Roll_Y - Roll_X)
    Area = 0.5*Area
    
    return Area
 
#############################################################################
###         Fits a function to find approx point
#############################################################################
def FitFct_Exp(XData, YData, UsePoints=None, plot=False, **kwargs):
#############################################################################  
    """
    Fits an higher order order polynom function against a set of data in regime 
    around a fitting point

    paramters              description
    =====================  =============================================:
    XData                   original X-Axis Data
    YData                   original Y-Axis Data
    plot                    (option) plot results/fitting for debug
    
    return type
       approximated a + b*np.exp(c*x) values
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        XData = ...
        YData = ...
        
        
        # Find Point
        [a,b] = pcb.FindPoint_FitFct(XData, YData,)
           
   
    """
#############################################################################    
    
    def parse_slice(value):
        """
        Parses a `slice()` from string, like `start:stop:step`.
        """
        if value:
            parts = value.split(':')
            if len(parts) == 1:
                # slice(stop)
                parts = [None, parts[0]]
            # else: slice(start, stop[, step])
        else:
            # slice()
            parts = []
        return slice(*[int(p) if p else None for p in parts])

#############################################################################   
    # Fit Function
    def Fitting_func(x, a, b, c):
        return a + b*np.exp(c*x)
    
    # cut data?
    if UsePoints:
        XData_Use = XData[parse_slice(UsePoints)]
        YData_Use = YData[parse_slice(UsePoints)]
    else:
        XData_Use = XData
        YData_Use = YData
        
    # fit curve
    popt, pcov = curve_fit(Fitting_func,
                           XData_Use, 
                           YData_Use,
                           **kwargs)
    
    # Error
    if np.any(pcov == np.inf):
        print("Fitting Error!")
    
    
    if plot==True:
        
        # plot fitted variables
        print('fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
        
        # fitted curve:
        Yfit = Fitting_func(XData, *popt)
        
        plt.figure()
        plt.plot(XData,YData, 'x')
        plt.plot(XData,Yfit, 'x')
        plt.ylim([0,200])
        plt.show
    
    return popt

#############################################################################
###         Fits a function to find approx point
#############################################################################
def FindPoint_FitFct(XData, YData, XPoint, Approx_Range, plot=False, **kwargs):
#############################################################################  
    """
    Fits an higher order order polynom function against a set of data in regime 
    around a fitting point

    paramters              description
    =====================  =============================================:
    XData                   original X-Axis Data
    YData                   original Y-Axis Data
    XPoint                  Fitting Point
    Approx_Range            range around point for fitting
    plot                    (option) plot results/fitting for debug
    
    return type
       approximated Y-Value for that point
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        XData = ...
        YData = ...
        
        XPoint = ...
        
        
        # Find Point
        YPoint = pcb.FindPoint_FitFct(XData, 
                                      YData, 
                                      XPoint, 5)
           
   
    """
#############################################################################    
    
    # Fit Function
    def Fitting_func(x, a, b, c, d):
        return a + b*x + c*x**2 + d*x**3
    
    # Find nearest point
    Xindex = np.argmin(np.abs(XData-XPoint))
    
    # approx range to large?
    if Approx_Range > Xindex:
        Approx_Range = Xindex
    
    # limit approx range
    Index_min = Xindex - Approx_Range
    Index_max = Xindex + Approx_Range
    
    # fit curve
    popt, pcov = curve_fit(Fitting_func,
                           XData[Index_min:Index_max],
                           YData[Index_min:Index_max],
                           **kwargs)
    
    # Error
    if np.any(pcov == np.inf):
        print("Fitting Error!")
    
    # Find Y Point
    YPoint = Fitting_func(XPoint, *popt)
    
    if plot==True:
        
        # plot fitted variables
        print('fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
        
        # fitted curve:
        Yfit = Fitting_func(XData, *popt)
        
        plt.figure()
        plt.plot(XData,YData, 'x')
        plt.plot(XData,Yfit, 'x')
        plt.show
    
    return YPoint

#############################################################################
###         Fits a function to find approx point
#############################################################################
def Linearization_Point(XData, YData, XPoint, Tolerance, num=100, 
                        iterationmax = 100, plot=True, **kwargs):
#############################################################################  
    """
    Fits an higher order polynom function against a set of data in regime 
    around a fitting point to linearize this point

    paramters              description
    =====================  =============================================:
    coordinate_file         csv file
    draw_polygon            (option) plot the input polygon
    delimiter               (option) csv file delimiter
    
    return type
       float with the calculated area
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # import csv
        file = 'test.csv'
        
        # calculate and print
        Area = pcb.CSV2Area(file, draw_polygon=True)
        
        print('This Polygon has an Area of ' + str(Area) + ' mm2')      
   
    """
#############################################################################    
    
    # Fit Function
    def FirstOrder_func(x, a, b):
        return a + b*x
    
    def ReferenceOrder_func(x, a, b, c, d, e, f, g):
        return a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5 + g*x**6
    
    # Find nearest point
    Xindex = np.argmin(np.abs(XData-XPoint))
    
    # Starting Point for index
    index_up = 1
    index_down = 1
    
    # Linearize around this point
    Index = [Xindex - index_up, Xindex + index_down]

    
    tolerance_borders = -1
    Index_Toleranceband = 0
    iteration = 0
 
############################################################################# 
    
    while (tolerance_borders%2 !=0):

        # limit iterations
        iteration = iteration+1
        if iteration == iterationmax:
            print("max iterations reached") 
            break
        
        # fit curve first order
        popt_lin, pcov = curve_fit(FirstOrder_func,
                               XData[Index[0]:Index[1]],
                               YData[Index[0]:Index[1]],
                               **kwargs)
        
        #fit curve third order against complete curve
        popt_3Order, pcov = curve_fit(ReferenceOrder_func,
                                      XData,
                                      YData,
                                      **kwargs)    
        
        # finer x_scale
        X_fine = np.linspace(XData[Index[0]], XData[Index[1]], num)
    
        # Calculate Error
        Yfit_1 = FirstOrder_func(X_fine, *popt_lin)
        Yfit_3 = ReferenceOrder_func(X_fine, *popt_3Order)
        Lin_Error = (Yfit_3-Yfit_1)/Yfit_3
        Lin_Error_Tolerance = np.abs(Lin_Error) - Tolerance
        Index_Toleranceband = np.where(np.diff(np.signbit(Lin_Error)))[0]

        # get iterations
        tolerance_borders = np.size(Index_Toleranceband)
        
        if plot==True:
            print("searching for tolerance borders. Currently found: "
                  + str(tolerance_borders) 
                  + " at Position "
                  + str(Index_Toleranceband)
                  + ". Current Limits: "
                  + str(Index))

        # Nothing Found?
        if tolerance_borders == 0:
            tolerance_borders = -1
            
            # Increase both indicies
            Index[0] = Index[0] - 1
            Index[1] = Index[1] + 1
            
            
        # Only one found?
        if tolerance_borders == 1:
            
            # Lower Border found? -> Increase upper one
            if Index_Toleranceband[0] < np.size(X_fine)/2:
                Index[1] = Index[1] + 1
                
            # Upper Border found? -> Decrease lower one
            if Index_Toleranceband[0] >= np.size(X_fine)/2:
                Index[0] = Index[0] - 1
                
         # Three found?
        if tolerance_borders == 3:
            
            # Last one Border found? -> Increase upper one
            if Index_Toleranceband[-1] < np.size(X_fine)*3/4:
                Index[1] = Index[1] + 1
                
            # Upper Border found? -> Decrease lower one
            if Index_Toleranceband[0] >= np.size(X_fine)/4:
                Index[0] = Index[0] - 1
                
        # Index in limits?
        if Index[0] < 0 or  Index[0] >= np.size(XData):
            print("No Convergence with parameters on this data set")
            iteration = iterationmax
            break
                   
        if Index[1] < 0 or Index[1] >= np.size(XData): 
            print("No Convergence with parameters on this data set")
            iteration = iterationmax
            break


#############################################################################               
    # Iterations done
    
    
    # use last values
    if iteration == iterationmax: 
        Index_Toleranceband = [0,np.size(X_fine)]
    
    # Fitted Dataset inside of Tolerance
    X_Lin = X_fine[Index_Toleranceband[0]:Index_Toleranceband[-1]]
    Y_Lin = FirstOrder_func(X_Lin, *popt_lin)   
    
    if plot==True:
        
        # plot fitted variables
        print('fit: a=%5.3f, b=%5.3f' % tuple(popt_lin))
        print(str(Index_Toleranceband))
          
        Yfit_3_complete = ReferenceOrder_func(XData, *popt_3Order)    
            
        plt.figure(figsize=(10,5))
        ax1 = plt.subplot(121)
        ax1.plot(XData,YData, 'x-', label="Data")
        ax1.plot(XData,Yfit_3_complete, '--', label="Reference Approx")   
        ax1.plot(X_fine,Yfit_1, ':', label="1st Order Approx", linewidth=4)   
        ax1.plot(X_Lin,Y_Lin, label="Linearized", linewidth=4)
        ax1.grid()
        ax1.legend()
        
        ax2 = plt.subplot(122)
        ax2.plot(X_fine,Lin_Error_Tolerance, label="rel. Error")   
        ax2.plot(X_fine,Lin_Error, label="|rel. Error| - Tolerance")           
        ax2.grid()
        ax2.legend()
        plt.show
    
    return [X_Lin, Y_Lin]

#############################################################################
###         Generate Plot for Time Domain / Linear
#############################################################################
def Linear_Plot(ax, Plot_list, X_label, Y_label, Legend=True, LegendLoc=0,
                TwinX=None, Ylim=None, XAutolim=True, fontsize=12, **kwargs):
#############################################################################  
    """
    Prepares a X-Y linear plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    Plot_list               all X and Y Values also Labels (and arguments)
    X_label                 X Axis Label and Unit (Engineering Package)
    Y_label                 Y Axis Label and Unit (Engineering Package)
    Legend                  (option) plot legend
    LegendLoc               (option) legend location
    TwinX                   (option) primary Y-Axis
    Ylim                    (option) set Y-Axis limits
    Xlim                    (option) set automatically X Limit
    fontsize                (option) Fontsize of this Plot 
    
    return type
       None  (writes directly into axis)
      
    Example:
        
        import PCB_Toolbox as pcb
        
        ...
        
        # Prepare
        Xlabel = [YLabel, 'V']
        Ylabel = [YLabel, 'V']
        Plot = [[XData, YData, "Label"], [XData, YData, "Label", 'linestyle=dashed'],...]   
        
        # Generate Plot
        plt.figure(figsize=(7.5,12))
        ax1 = plt.subplot(111)
        pcb.Linear_Plot(ax1, Plot, X_label, Ylabel)  
        plt.show()
   
    """        
#############################################################################      
    for plot in Plot_list:
        
        # check dimension of X-Axis if whole trace
        x_plot = plot[0]
        
        # only one marker?
        if np.size(x_plot) > 1:
            y_plot = plot[1][0:np.size(x_plot)]
        else:
            y_plot = plot[1]
        
        # emtpy argument list
        userargs = {}
        
        # insert plotting arguments
        if len(plot) >= 4:
            userargs = dict(e.split('=') for e in plot[3].split(', '))
            
        ax.plot(x_plot, y_plot, label=plot[2], **userargs)

    # label
    ax.set_ylabel(Y_label[0])
    ax.yaxis.set_major_formatter(tck.EngFormatter(unit=Y_label[1]))
    ax.set_xlabel(X_label[0])
    ax.xaxis.set_major_formatter(tck.EngFormatter(unit=X_label[1]))
    
    # xlimit
    if XAutolim:
        
        # search min and max x values
        x_limit_min = np.min(plot[0])
        x_limit_max = np.max(plot[0])
        
        # iterate all traces
        for trace in ax.get_lines():
            
            # find new min value
            if np.min(trace.get_xdata()) < x_limit_min:
                x_limit_min = np.min(trace.get_xdata())
 
            # find new min value
            if np.max(trace.get_xdata()) > x_limit_max:
                x_limit_max = np.max(trace.get_xdata())           
        
        # set x limit
        ax.set_xlim([x_limit_min,x_limit_max])
    
    # xlimit    
    if Ylim:
        ax.set_ylim([Ylim[0],Ylim[1]])
        
    # set font sizes
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
    
        # Align both Y Axis to grid
    if not(TwinX==None) and type(TwinX) == type(ax) and Legend:
        
        # Align Axis
        Align_YAxis(ax,TwinX)
        
        # include axis labels in single legend
        all_lines = TwinX.get_lines() + ax.get_lines()
        all_labels = [l.get_label() for l in all_lines]
        
        # plot legend
        TwinX.legend(all_lines, all_labels, framealpha=1, loc=LegendLoc)
        
    else:
        
        if Legend:
            # legend
            ax.legend(framealpha=1, loc=LegendLoc)
    
        # Generate new Grid
        ax.minorticks_on()
        ax.grid(True, which='major')
        ax.grid(which='minor', alpha=1, linestyle=':', linewidth=1)
        ax.grid(which='major', alpha=1, linewidth=1.2)    
        
    #retrn
    return ax

#############################################################################
###         Generate Plot for Frequency Domain / SemilogX
#############################################################################
def SemiLogX_Plot(ax, Plot_list, X_label, Y_label, Legend=True, LegendLoc=0, 
                  TwinX=None, Ylim=None, XAutolim=True, fontsize=12, 
                  LegendAlpha=1, **kwargs):
#############################################################################  
    """
    Prepares a X Log and Y Linear plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    Plot_list               all X and Y Values also Labels (and arguments)
    X_label                 X Axis Label and Unit (Engineering Package)
    Y_label                 Y Axis Label and Unit (Engineering Package)
    Legend                  (option) plot legend
    LegendLoc               (option) legend location
    TwinX                   (option) primary Y-Axis
    Ylim                    (option) set Y-Axis limits
    Xlim                    (option) set automatically X Limit
    fontsize                (option) Fontsize of this Plot 
    LegendAlpha             (option) Transparency of legend box
    
    return type
       None  (writes directly into axis)
       
    Example:
        
        import PCB_Toolbox as pcb
        
        ...
        
        # Prepare
        Xlabel = [YLabel, 'V']
        Ylabel = [YLabel, 'V']
        Plot = [[XData, YData, "Label",, 'linestyle=dashed']]    
        
        # Generate Plot
        plt.figure(figsize=(7.5,12))
        ax1 = plt.subplot(111)
        pcb.Linear_Plot(ax1, Plot, X_label, Ylabel)  
        plt.show()
   
    """    

#############################################################################      
 
    for plot in Plot_list:
        
        # check dimension of X-Axis
        x_plot = plot[0]
        
        # only one marker?
        if np.size(x_plot) > 1:
            y_plot = plot[1][0:np.size(x_plot)]
        else:
            y_plot = plot[1]
        
        # emtpy argument list
        userargs = {}
        
        # insert plotting arguments
        if len(plot) >= 4:
            userargs = dict(e.split('=') for e in plot[3].split(', '))

        # plot
        ax.semilogx(x_plot, y_plot, label=plot[2], **userargs)     
       
    # label
    ax.set_ylabel(Y_label[0])
    ax.yaxis.set_major_formatter(tck.EngFormatter(unit=Y_label[1]))
    ax.set_xlabel(X_label[0])
    ax.xaxis.set_major_formatter(tck.EngFormatter(unit=X_label[1]))
        
    # set font sizes
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
        
        # Align both Y Axis to grid
    if not(TwinX==None) and type(TwinX) == type(ax) and Legend:
        
        # Align Axis
        Align_YAxis(ax,TwinX)
        
        # include axis labels in single legend
        all_lines = TwinX.get_lines() + ax.get_lines()
        all_labels = [l.get_label() for l in all_lines]
        
        # plot legend
        TwinX.legend(all_lines, all_labels, framealpha=LegendAlpha, loc=LegendLoc)
        
    else:
        
        if Legend:
            # legend
            ax.legend(framealpha=1, loc=LegendLoc)
    
        # Generate new Grid
        ax.minorticks_on()
        ax.grid(True, which='major')
        ax.grid(which='minor', alpha=1, linestyle=':', linewidth=1)
        ax.grid(which='major', alpha=1, linewidth=1.2)
        
        # xlimit
    if XAutolim:
        ax.set_xlim([np.min(plot[0]),np.max(plot[0])])

    # xlimit    
    if Ylim:
        ax.set_ylim([Ylim[0],Ylim[1]])
        
    # jump back
    return

#############################################################################
###         Generate Vertical Line with Label
#############################################################################
def Vline_Plot(ax, xValue, xLabel, yDistance=0.25, yPos='up', color='r',
               fontsize='12', linestyle='-', horizontalalignment='center',
               **kwargs):
#############################################################################  
    """
    Generates Vertical Line in Plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    xValue                  Value on X-Axis
    xLabel                  Label for Line
    yDistance               Distance Factor (Y Tick Distance, default=0.25)
    yPos                    'up' or 'down'
    color                   color of line and text (default=red)
    fontsize                fontsize of text (default=12)
    linestyle               linestyle of line (default='-')
    horizontalalignment     Alignment of text (default='center')
    
    return type
       None  (writes directly into axis)
       
    Example:
        
        import PCB_Toolbox as pcb

   
    """    

#############################################################################  
    # Add vertical line
    ax.axvline(x=xValue, color=color, linestyle=linestyle)

    # find y Position
    ylimits = ax.get_ylim()
    ydistance = np.mean(np.ediff1d(ax.get_yticks()))*yDistance    
    
    # up or down?
    if yPos == 'up':
        ylimits = ylimits[1]
        ydistance = np.abs(ydistance)
        
    if yPos == 'down':
        ylimits = ylimits[0]
        ydistance = -1*np.abs(ydistance) 

    # generate Text            
    ax.text(xValue, ylimits+ydistance, xLabel, color=color,
            fontsize=fontsize, horizontalalignment=horizontalalignment)  
    
    # jump back
    return
    
#############################################################################
###         Generate Vertical Line with Label
#############################################################################
def Hline_Plot(ax, yValue, yLabel, xDistance=0.4, xPos='right', color='r',
               fontsize='12', linestyle='-', verticalalignment='center',
               **kwargs):
#############################################################################  
    """
    Generates Vertical Line in Plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    yValue                  Value on y-Axis
    yLabel                  Label for Line
    xDistance               Distance Factor (x Tick Distance, default=0.25)
    xPos                    'left' or 'right'
    color                   color of line and text (default=red)
    fontsize                fontsize of text (default=12)
    linestyle               linestyle of line (default='-')
    verticalalignment       Alignment of text (default='center')
    
    return type
       None  (writes directly into axis)
       
    Example:
        
        import PCB_Toolbox as pcb

   
    """    

#############################################################################  
    # Add vertical line
    ax.axhline(y=yValue, color=color, linestyle=linestyle)

    # find y Position
    xlimits = ax.get_xlim()
    xdistance = np.mean(np.ediff1d(ax.get_xticks()))*xDistance    
    
    # up or down?
    if xPos == 'right':
        xlimits = xlimits[1]
        xdistance = np.abs(xdistance)
        
    if xPos == 'left':
        xlimits = xlimits[0]
        xdistance = -1*np.abs(xdistance) 

    # generate Text            
    ax.text(xlimits+xdistance, yValue, yLabel, color=color,
            fontsize=fontsize, verticalalignment=verticalalignment)  
    
    # jump back
    return
     
#############################################################################
###         Align two Y-axis
#############################################################################
def Align_YAxis(ax1, ax2):
#############################################################################    
    """
   Align two Axis on the same Grid

    paramters              description
    =====================  =============================================:
    ax1,ax2                 axis object
    
    return type
       two vectors for axis ticks
    
    """
#############################################################################     

    # Generate new Ticks Function
    def Generate_newTicks(ax_dy, ax_ystart, max_ticks):
        
        # Check Intervall
        ax_intervall = np.ceil(ax_dy/(max_ticks - 1))
        
        # Generate new Intervall
        ax_dy_new = (max_ticks -1) * ax_intervall
        
        # generate new Ticks                
        ax_ticks_new = np.linspace(ax_ystart,
                                   ax_ystart+ax_dy_new,
                                   max_ticks)
        
        return ax_ticks_new

#############################################################################  
        
    # Find Round Digits
    def Generate_RoundPoint(ax_dy, max_ticks):
        
        # Check Decimal
        log = np.log10(abs(ax_dy/max_ticks))
        exponent = np.floor(log)-1
        
        # Generate Rounded Value
        ax_roundto = 10**exponent
        
        return ax_roundto
    
#############################################################################  

    # get maximum number of ticks
    max_ticks = max(len(ax1.get_yticks()),len(ax2.get_yticks()))
            
    # get axis scaling
    ax1_dy = ax1.get_ybound()[1] - ax1.get_ybound()[0]   
    ax2_dy = ax2.get_ybound()[1] - ax2.get_ybound()[0]  

    # Roundto Number
    ax1_roundto = Generate_RoundPoint(ax1_dy, max_ticks)
    ax2_roundto = Generate_RoundPoint(ax2_dy, max_ticks)
     
    # get axis bounds and scale
    YBound_ax1 = [np.floor(ax1.get_ybound()[0]/ax1_roundto),
                  np.ceil(ax1.get_ybound()[1]/ax1_roundto)]
    YBound_ax2 = [np.floor(ax2.get_ybound()[0]/ax2_roundto),
                  np.ceil(ax2.get_ybound()[1]/ax2_roundto)]

    # get axis scaling and scale
    ax1_dy = YBound_ax1[1] - YBound_ax1[0]   
    ax2_dy = YBound_ax2[1] - YBound_ax2[0]
    
    # define starting points
    ax1_start = YBound_ax1[0]
    ax2_start = YBound_ax2[0]
     
    # Generate new Ticks for axis
    ax1_ticks_new =   Generate_newTicks(ax1_dy, ax1_start, max_ticks)
    ax2_ticks_new =   Generate_newTicks(ax2_dy, ax2_start, max_ticks)
    
    # reverse rescaling
    ax1_ticks_new = ax1_ticks_new*ax1_roundto
    ax2_ticks_new = ax2_ticks_new*ax2_roundto
    
    # set ticks
    ax1.set_yticks(ax1_ticks_new)
    ax2.set_yticks(ax2_ticks_new)
    
  
    # jump back    
    return

#############################################################################
###         Fits a function to find approx point
#############################################################################
def FindPoint_NextValue(XData, YData, XPoint, plot=False, **kwargs):

#############################################################################  
    """
    Fits an higher order order polynom function against a set of data in regime 
    around a fitting point

    paramters              description
    =====================  =============================================:
    XData                   original X-Axis Data
    YData                   original Y-Axis Data
    XPoint                  Fitting Point
    
    return type
       give nearest minimum value
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        XData = ...
        YData = ...
        
        XPoint = ...
        
        
        # Find Point
        YPoint = pcb.FindPoint_NextValue(XData, YData, XPoint)
           
    """
#############################################################################     
    # Search next minimum
    index = np.argmin(np.abs(XData-XPoint))

    return YData[index]

#############################################################################
###         Digitalize Data
#############################################################################
def Digitalize_Data(data, clock, edge_trigger = 'rising',
                    high_val = 1, low_val = 0, trigger_val = 0.5,
                    threshold_high = 2, threshold_low = 0.8, plot = False):
#############################################################################  
    """
    Fits an higher order order polynom function against a set of data in regime 
    around a fitting point

    paramters              description
    =====================  =============================================:
    data                    data array
    clock                   clock array
    edge_trigger            (optional) rising or falling
    high_val                (optional) parsed high value
    low_val                 (optional) parsed low value
    trigger_val             (optional) trigger value
    threshold_high          (optional) threshold high value
    threshold_low           (optional) threshold low value
    plot                    (optional) plot data for debug
    
    return type
       binary bit stream
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        data = ...
        clock = ...
                
        
        # Find Point
        stream = pcb.Digitalize_Data(data, clock)
           
    """   
#############################################################################   
     
    # Schmitt-Trigger function
    def digitialize(analog_array, high_val, low_val, threshold_high, threshold_low):
        
        # empty variables
        conversion_errors = []
        digital_array = []
    
        # digitalize data in array
        current_index = 0
        for analog_value in analog_array:
            
            # Schmitt trigger
            if analog_value >  threshold_high:
                digital_array.append(high_val)
                
            elif analog_value < threshold_low:
                digital_array.append(low_val)
                
            # out of range
            else:
               conversion_errors.append(current_index)
               digital_array.append(analog_value)
               print("Data " + str(analog_value) + " (index = " 
                     + str(current_index) + ") could not be parsed!") 
                
            # increase index
            current_index = current_index + 1
            
        return [digital_array, conversion_errors]
    
#############################################################################    
            
    # Generate digital stream
    data = digitialize(data, high_val=high_val, low_val=low_val,
                               threshold_high=threshold_high, threshold_low=threshold_low)
    
    clock = digitialize(clock, high_val=high_val, low_val=low_val,
                               threshold_high=threshold_high, threshold_low=threshold_low)  
    
    # Remove Indicies
    conversion_errors = data[1] + clock[1]
    
    # clear faulty data from both streams
    data = [value for index, value in enumerate(data[0]) if index not in conversion_errors]
    clock = [value for index, value in enumerate(clock[0]) if index not in conversion_errors]    
 
    # data conversion
    data = np.array(data)
    clock = np.array(clock)
     
    # rising or falling edge clock mask
    mask_rising = (clock[:-1] < trigger_val) & (clock[1:] > trigger_val)
    mask_falling = (clock[:-1] > trigger_val) & (clock[1:] < trigger_val) 
    
    # choose edge mask
    mask = []
    if edge_trigger == 'rising':
        mask = [index for index, value in enumerate(mask_rising) if value]
    elif edge_trigger == 'falling':
        mask = [index for index, value in enumerate(mask_falling) if value]
    else :
        print("edge trigger mask not valid!")         

    # sample data at edge
    binary_data = data[mask]
           
    # plot both streams
    if plot:
             
        plt.figure(figsize=(15,10))
        ax1 = plt.subplot(311)
        ax1.plot(data, label="Data")
        ax1.plot(mask,binary_data, 'rx', label="Sampled")
        ax1.grid()
        ax1.legend()
        
        ax2 = plt.subplot(312)
        ax2.plot(clock, label="Clock")           
        ax2.grid()
        ax2.legend()
        
        ax3 = plt.subplot(313)
        ax3.plot(mask_rising, label="Rising Edge")  
        ax3.plot(mask_falling, label="Falling Edge")           
        ax3.grid()
        ax3.legend()
        plt.show       
             
    return binary_data
 
#############################################################################
##          Convert Complex Data into Polar
#############################################################################
def CMPLX2Format (complex_data, voltage=False):
#############################################################################  
    """
    Converts complex data (a+jb) into different formats (polar, dB, ..)

    paramters              description
    =====================  =============================================:
    complex_data            input data
    voltage                 (optional) convert in dB with 20*log10(..)
    
    return type
       struct with complex, dB, mag, phase(degree)
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        data = ...               
        
        # Generate multiple formats 
        data = pcb.CMPLX2Polar(data)
           
    """   
#############################################################################  
    # generate output struct
    export_data = {}
       
    # complex data
    export_data['complex'] = complex_data
    
    # convert Magniture in dB
    if voltage:
        export_data['dB'] = 20*np.log10(np.abs(complex_data))
    else:
        export_data['dB'] = 10*np.log10(np.abs(complex_data))
        
    # convert Magniture linear
    export_data['mag'] = np.abs(complex_data)
    
    # convert Phase in Degree
    export_data['PhaseDeg'] = np.unwrap(np.angle(complex_data, deg=True))
    
    return export_data

#############################################################################
##          Time Domain Averaging
#############################################################################
def Average(XData, YData, Points=200):
############################################################################# 
    """
    Average of YData, both data arrays will be splitted into eqal chunks

    paramters              description
    =====================  =============================================:
    YData                   input data, which should be averaged
    XData                   input data, average data will correspondent to middle
    Points                 (optional) number of points per chunk
    
    return type
       struct with List YData and XData
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        dataY = ...               
        dataX = ...       
          
        # Generate multiple formats 
        dataAVG = pcb.Average(dataY, dataX)
           
    """   
############################################################################# 
    
    # calculate chunk section number
    Chunk_AllIndices = np.int(np.size(YData)/Points)

    # Split Array
    YData = np.array_split(YData, Chunk_AllIndices)
    XData = np.array_split(XData, Chunk_AllIndices)

    # Filtering
    for Chunk_index in np.arange(Chunk_AllIndices):
        
        # Average of Dataset
        YData[Chunk_index] = np.average(YData[Chunk_index])
        
        # Generate new X Row
        XData[Chunk_index] = XData[Chunk_index][0]+((XData[Chunk_index][-1] - XData[Chunk_index][0])/2)

    # Return Data  
    return [XData, YData]

#############################################################################
##          FFT of Data Signal
#############################################################################
def FFTData(XData, YData, ZeroPad=1000):
############################################################################# 
    """
    FFT of input data
    
    paramters              description
    =====================  =============================================:
    YData                   input data, which should be averaged
    XData                   input data, average data will correspondent to middle
    ZeroPad                 (optional) zeropadding with numbers
    
    return type
       struct with List YData and XData
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        dataY = ...               
        dataX = ...       
          
        # Generate multiple formats 
        dataAVG = pcb.Average(dataY, dataX)
           
    """   
############################################################################# 
    # ZeroPadding
    YData = np.pad(YData, ZeroPad, mode='constant')
    
    # estimate Sample Time (mean)
    Tsample = np.mean(XData[1:-1] - XData[0:-2])

    # Number of Samples
    NFFT = np.size(YData)

    # generate FFT and power spectrum 
    yFFT= np.fft.fft(YData)
    
    # Rescale Amplitude and only positive frequency axis
    yFFT = 2.0/NFFT * np.abs(yFFT[:(NFFT//2)])

    # Generate Frequency Axis
    xFFT = np.linspace(0,1/(2*Tsample),NFFT/2)
    
    return [xFFT, yFFT]

#############################################################################
##          FFT of Data Signal
#############################################################################
def String2List(file, delimiter=','):
############################################################################# 
    """
    Average of YData, both data arrays will be splitted into eqal chunks

    paramters              description
    =====================  =============================================:
    file                    input data, which should be averaged
    
    return type
       list with data
       
    Example:
        
        import PCB_Toolbox as pcb
        
        # generate Dataset
        file = ...  
             
        # Generate multiple formats 
        list = pcb.String2List(file)
           
    """   
############################################################################# 
    # Generate a List with strings of data
    Input_lines = [line.rstrip('\n') for line in open(file)]

    Output = []

    # Generate of List in List
    for line in Input_lines:
        
        try:
           # Seperate at delimiter
            Output.append([float(s) for s in line.split(delimiter)])
                    
        except:
            
            # could not convert
            Output.append(line)


    return Output




    