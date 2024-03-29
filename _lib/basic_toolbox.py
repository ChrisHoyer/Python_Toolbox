#############################################################################
#
# - CSV2Dict: script for importing CSV Files
# - CSV2Area: script for calculating polygon area based on coordinates
# - FitFct_Exp: script that fits data against an Exp Function
# - FindPoint_FitFct: fits higher order Polynom against Dataset to find point
# - Linearization_Point: linearization around one point
# - Linear_Plot: linear plot function with automated labeling
# - Polar_Plot: Polar plot function with automated labeling
# - Histogram_Plot: Plot Histogram
# - Box_Plot: Boxplots with automated labeling
# - SemiLogX_Plot: semilog x plot function with automated labeling
# - Vline_Plot: generates vertical line in plot with label
# - Hline_Plot: generates horizontal line in plot with label
# - Fill_Plot: generates fill pattern between axis
# - Rectangle_Plot: generates rectangle inside plot
# - Align_YAxis: Align two YGrids in one plot
# - FindPoint_NextValue: find nearest point
# - Digitalize_Data: Generates binary stream from data and clock from
# - Extract_DigitalClock: Extracts clock from data and clock from waveform
# - CMPLX2Format: converts complex numbers in different formats
# - Average: Average a Dataset
# - FrequencyFiltering: Frequency Filtering of Data
# - String2List: Generates a List out of a csv with one line
# - XYZ_Plot: generates 3D Plot (X,Y,Z) for i.e. waterfall diagrams
# - Spectrum_Minimizer: Generates Mean and Max Spectrum from Dataset
# - MovingFilter: Moving Average, Max, Min Filter
# - printProgressBar: Print progress Bar in Console
# - EngNot: Engineering Notation (Si-Prefix)
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Date: 06-11-2021
#############################################################################

from scipy.optimize import curve_fit
from cycler import cycler

import matplotlib.ticker as tck
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import csv
import re


#############################################################################
###         Global Variables
#############################################################################

# Black and White Style
monochrome = (cycler('color', ['k']) * cycler('linestyle', ['-', '--', ':']) * cycler('marker', ['^', '.','v', '<', '>']))


# Grid Crosses
# https://stackoverflow.com/questions/54132998/draw-markers-at-intersections-of-grids-minor-ticks
def set_grid_cross(ax):
    xticks = ax.get_xticks(minor=True)
    yticks = ax.get_yticks(minor=True)
    xgrid, ygrid = np.meshgrid(xticks, yticks)
    kywds = dict() 
    grid_lines = ax.plot(xgrid, ygrid, 'o', markersize=0.1, color='lightgray', alpha=0.5)
    
#############################################################################
###         Import CSV File to a dictionary
#############################################################################
def CSV2Dict(csv_file, delimiter=';', complexdelimiter='/',  Multiblock = False,
             headerline_ends=0, blockheader=False, blockkeys=[],
             headerkeys=[], cellsize=[], **kwargs):
############################################################################# 
    """
    Imports all CSV Data (also complex and logarithmic data)
    Specialized for LTspice and ADS

    paramters              description
    =====================  =============================================:
    csv_file                csv file
    delimiter               (option) csv file delimiter
    complexdelimiter        (option) delimiter if there is a complex value
    Multiblock              (option) There are multiple lists nested in one file
    blockheader             (option) each block has header?
    headerline_ends         (option) end of header (skip until this line, -1: skip header at all)
    blockkeys               (option) custom block header names (must be equal to block count)
    headerkeys              (option) custom header names (must be equal to row count)
    cellsize                (option) cellsize of each cell
    
    return type
       dictionary with each column or nested dictionary
       
    Example:
        
        import basic_toolbox as basic
        
        # import csv
        file = 'test.csv'
        
        # calculate
        Data = basic.CSV2Dict(file)
         
   
    """
#############################################################################  
    # row state
    States = ['', 'SimpleNumber', 'ComplexNumber', 'Text', 'Array', 'EmptyRow']
    
    # Matching Numbers to parse float / scientific data
    scientific_number = re.compile('[-+]?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?')  
    
            
#############################################################################
    
    def CSVData_DefineContent(header, row, v):
        
        # Default: No State
        Current_State = States[0]
        
        # remove leading or tailing whitespaces
        cell_content = v.strip()
                
        # only plain Number (normal or scientific)
        try:
            float(cell_content)
            Current_State = States[1]
                
        except:
             
            # Row with Text and empty
            if (type(cell_content) is str) & (cell_content == "") :
                 Current_State = States[5]  
                 
            # Row with Text and not empty?
            elif (type(cell_content) is str) & (cell_content != "") :
                 Current_State = States[3]  
                 
            # Data is Array?
            elif (type(cell_content) is str) & (cell_content[0] == '[') & (cell_content[-1] == ']'):
                 Current_State = States[4]                     
 
            # Contains two values? -> Complex
            elif len(cell_content.split(complexdelimiter)) == 2:
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

            # Convert Text Data                                
            elif State == States[3]:
                v_parsed = v
                            
            # Convert Header Data                                
            elif State == States[5]:
                # new Block starts
                v_parsed = States[5]
                
            # Array?
            elif State == States[4]:
                
                # generate list
                v_parsed = v.replace('[', "").replace(']', "")
                v_parsed = v_parsed.split(",")
                
                # float
                v_parsed = [float(i.strip()) for i in v_parsed]
                
            else:
                v_parsed = float('NaN')
                print("Data " + str(v) + " could not be parsed!")                  
            
        except:
            print("Data " + str(v) + " could not be parsed!") 
            
        return v_parsed
        
#############################################################################                

    # Import CSV Matrix
    csv_reader = csv.reader(open(csv_file), delimiter=delimiter)
		
		# start with first line
    global_line_count = 0
    
    # generate new Block
    NewBlockStart = True
    block_line = 0
    block_name = ''
    block_count = 0
    
    # Content
    Block = {}
    Data = {}
    
    
 #############################################################################  
    
    # Iterate each row of the file
    for row in csv_reader:
         
        # Skip Global Header
        if  global_line_count < (headerline_ends):
            global_line_count += 1
            continue
			 
        # empty row? -> skip and  new Block
        if not row:
            NewBlockStart = True
            continue
        
        # name new block
        if NewBlockStart:
            
            # reset flag
            NewBlockStart = False
            
            # Any external names defined?
            if not blockkeys:
                block_name = block_count
            else:
                block_name = blockkeys[block_count]
                
            # Insert readed Block and start new block
            if len(Block) > 0:
                
                # insert into dataset
                Data[block_name] = Block
                block_count = block_count + 1
                
                # delete old Block
                block_line = 0
                Block = {}
        
        # starting new block
        if block_line == 0:
            
            # Using current line as Header
            header = [r.strip() for r in row]
            
            # TODO? What happend, when there is no header? i.e. ADS export
            
            # external header defined and same size?
            if len(header) == len(headerkeys):
                header = headerkeys
                
            # generate header for empty block dictinary
            for h in header:
                Block[h] = []
 
        # skip headerline, if headerline_ends == -1
        if headerline_ends == -1:
                block_line = 1
 
 ############################################################################# 
                
        # data content begins
        if block_line > 0:
            
            # row is larger then header?
            if len(row) > len(header):
                
                #rearrange row to match cellsize
                if len(cellsize) == len(header):
                    
                    # generate new row
                    new_row = []
                    lastsize = 0
                    
                    # iterate cellsize and rearrange
                    for size in cellsize:

                        # only one item?
                        if size == 1:
                            newcell =  row[lastsize]
                            newcell = str(newcell).replace("'","")
                            lastsize = size
                            
                        else:
                            # bring into new row
                            newcell = row[lastsize:(size+lastsize)]
                            newcell = str(newcell).replace("'","")
                            lastsize = size + 1
                            

                        
                        # include new cell into new row
                        new_row.append(newcell)
                        
                    # generate new row, with right length
                    row = new_row
                    
                else:
                    print("Please include cellsize!")
                    return 0
                
            # iterate row with repect to header
            for header_cell, cell in zip(header,row):
                
                # classify data and parse content
                parsed_Data = CSVData_ParseContent(header, row, cell)
                
                # Header found?
                if (parsed_Data == States[5]) and Multiblock:
                    NewBlockStart = True
                    
                    # TODO New Block
                    continue
                    
                # add Data into Block
                Block[header_cell].append(parsed_Data)
                   
        # next line
        block_line += 1
        global_line_count += 1
        
  #############################################################################                    
        
    # only one Block? 
    if block_count == 0:
        return Block
    else:
        
        # add last block, if not empty
        if len(Block) > 0:
            Data[block_name] = Block
            
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
        
        import basic_toolbox as basic
        
        # import csv
        file = 'test.csv'
        
        # calculate and print
        Area = basic.CSV2Area(file, draw_polygon=True)
        
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
        
        import basic_toolbox as basic
        
        # generate Dataset
        XData = ...
        YData = ...
        
        
        # Find Point
        [a,b] = basic.FindPoint_FitFct(XData, YData,)
           
   
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
        
        import basic_toolbox as basic
        
        # generate Dataset
        XData = ...
        YData = ...
        
        XPoint = ...
        
        
        # Find Point
        YPoint = basic.FindPoint_FitFct(XData, 
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
        
        import basic_toolbox as basic
        
        # import csv
        file = 'test.csv'
        
        # calculate and print
        Area = basic.CSV2Area(file, draw_polygon=True)
        
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
###         Generate Generic Plot ( Check Plot Calls!)
#############################################################################
def Generic_Plot(func, ax, Plot_list, X_label, Y_label, Legend=True, LegendLoc=0,
                TwinX=None, TwinY=None, TwinReuseTicks="BOTH",  Ylim=None, Xlim=None,
                XAutolim=True, fontsize=7, TicksEng=True, XTicksLabel=None,
                YTicksLabel=None,legendcol=1,fontsize_label=8, yaxis_pad=0, xaxis_pad=0, 
                BlackWhite=False, grid = True, minorgridalpha=.3, majorgridalpha=.6,
                legendalpha=1, funcReturn=False, grid_lw_major=1.2, grid_lw_minor=1.0,
                grid_zorder=-10,
                **kwargs):
#############################################################################  
    """
    Prepares a X-Y linear plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    Plot_list               all X and Y Values also Labels (and matplotlib arguments)
    X_label                 X Axis Label and Unit (option: rescaling factor) (Engineering Package)
    Y_label                 Y Axis Label and Unit (option: rescaling factor)(Engineering Package)
    Legend                  (option) plot legend
    LegendLoc               (option) legend location
    TwinX                   (option) secondary Y-Axis
    TwinY                   (option) secondary X-Axis
    TwinReuseTicks          (option) Methode of regenerating ticks (NONE | AX1 = use AX1 | BEST = fith both)
    Ylim                    (option) set Y-Axis limits [Y0,Y1]
    Xlim                    (option) set X-Axis limits [X0,X1]
    XAutolim                (option) set automatically X Limit (bool)
    TicksEng                (option) Enable Engineering Ticks
    XTicksLabel             (option) Label only ever nth tick on X
    YTicksLabel             (option) Label only ever nth tick on Y
    fontsize                (option) Fontsize of the legend and ticks
    fontsize_label          (option) Fontsize of the axis labels
    legendcol               (option) Legend Columns
    yaxis_pad               (option) move label to y-axis (padding)
    xaxis_pad               (option) move label to x-axis (padding)    
    BlackWhite              (option) Use Black and White Preset
    grid                    (option) Use grid
    
    return type
       None  (writes directly into axis)
      
    Example:
        
        import basic_toolbox as basic
        
        ...
        
        # Prepare
        Xlabel = ["XLabel", 'V']
        Ylabel = ["YLabel", 'V']
        Plot = [[XData, YData, "Label"], [XData2, YData2, "Label", 'linestyle=dashed'],
                [XData2, YData2, "Label", 'linestyle=dashed, color=k'], ...]
        
        # or for loops just use
        Plot.append([XData3, YData3, "Label", "**kwargs"])
        
        # Generate Plot
        plt.figure(figsize=(7.5,12))
        ax1 = plt.subplot(111)
        basic.Linear_Plot(ax1, Plot, Xlabel, Ylabel)  
        plt.show()
        
        # Two Y Axis one plot (similar to X-Axis)
        ax2 = ax1.twinx()
        basic.Linear_Plot(ax2, Plot, Xlabel, Ylabel, TwinX=ax1)         
   
    """        
#############################################################################   

    # BlackWhite Default Settings
    if BlackWhite:
        ax.set_prop_cycle(monochrome)
        
    # for multiple returns
    returnvals = []
        
    # check if Plot has entries
    if len(Plot_list) == 0:
        print("Nothing to Plot!")
        return
    
    for index in range(len(Plot_list)):
        
        plot = Plot_list[index]
        
        # Call Specific Plotting Function
        if funcReturn:
            ax, x_plot, returnval = func(ax, plot, Y_label, X_label)
            returnvals.append(returnval)
            
        else:
            ax, x_plot = func(ax, plot, Y_label, X_label)
        
    # label
    ax.set_ylabel(Y_label[0], labelpad=yaxis_pad)
    ax.set_xlabel(X_label[0], labelpad=xaxis_pad)
    
    # ticks in engineering formatter
    if TicksEng:
        ax.yaxis.set_major_formatter(tck.EngFormatter(unit=Y_label[1]))
        ax.xaxis.set_major_formatter(tck.EngFormatter(unit=X_label[1]))
    
    # xlimit
    if XAutolim:
        
        # search min and max x values
        x_limit_min = np.min(x_plot)
        x_limit_max = np.max(x_plot)
                       
        # iterate all traces
        for trace in ax.get_lines():
            
            # Length should be more than 2 points
            if len(trace.get_xdata()) > 2:
            
                # find new min value
                if np.min(trace.get_xdata()) < x_limit_min:
                    x_limit_min = np.min(trace.get_xdata())
     
                # find new min value
                if np.max(trace.get_xdata()) > x_limit_max:
                    x_limit_max = np.max(trace.get_xdata())   
        
        # set x limit
        ax.set_xlim([x_limit_min,x_limit_max])

    # xlimit    
    if Xlim:
        ax.set_xlim([Xlim[0],Xlim[1]])
        
    # ylimit    
    if Ylim:
        ax.set_ylim([Ylim[0],Ylim[1]])
        
    # set font sizes (all)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
    
    # set font size label
    for item in ([ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(fontsize_label)

    # =================================== 
    # change XTick Label Position
    if XTicksLabel:
        
        # change visibility of each Nth tick
        for (index,label) in enumerate(ax.xaxis.get_ticklabels()):
            if index % XTicksLabel != 0:
                label.set_visible(False)

    # change YTick Label Position
    if YTicksLabel:
        
        # change visibility of each Nth tick
        for (index,label) in enumerate(ax.yaxis.get_ticklabels()):
            if index % XTicksLabel != 0:
                label.set_visible(False)                
                
    # ===================================    
    # Legend and grid for two axis
    if not(TwinX==None) and type(TwinX) == type(ax):

        # include axis labels in single legend
        all_lines = TwinX.get_lines() + ax.get_lines()
        all_labels = [l.get_label() for l in all_lines]
        
        if Legend:        
            # plot legend
            TwinX.legend(all_lines, all_labels, framealpha=1, loc=LegendLoc) 
        
        # Align Axis
        Align_YXAxis(ax, TwinX, AxisType="Y", Method=TwinReuseTicks)
        
    elif not(TwinY==None) and type(TwinY) == type(ax):
 
        # include axis labels in single legend
        all_lines = TwinY.get_lines() + ax.get_lines()
        all_labels = [l.get_label() for l in all_lines]

        if Legend:        
            # plot legend
            TwinY.legend(all_lines, all_labels, framealpha=legendalpha, loc=LegendLoc) 
        
        # Align Axis
        Align_YXAxis(ax, TwinY, AxisType="X", Method=TwinReuseTicks)


    # ===================================    
    # grid and legend
    else:
        
        if Legend:
            # legend
            ax.legend(framealpha=legendalpha, loc=LegendLoc, fontsize=fontsize, ncol=legendcol)
            
        if grid:
            
            ax.minorticks_on()
            ax.grid(which='major', alpha=majorgridalpha, linestyle='-',
                    linewidth=grid_lw_major, zorder=grid_zorder) 
            ax.grid(which='minor', alpha=minorgridalpha, linestyle=':',
                    linewidth=grid_lw_minor, zorder=grid_zorder)
            
            

    #return
    if funcReturn:
        return returnvals
    else:
        return ax

#############################################################################
###         Generic Plotting Calls
#############################################################################

# pcolor Plot
def PseudoColor_Plot(ax, Plot_list, X_label, Y_label, **kwargs):

    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[XData, YData, ZData, "Label"],
                                 [XData2, YData2, ZData2, "Label", 'linestyle=dashed'],
                                 [XData2, YData2, ZData2, "Label", 'linestyle=dashed, color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """        
                
        # check dimension
        x_plot = plot[0]
        y_plot = plot[1]
        z_plot = plot[2]
        
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[4].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:

            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue
                
        # rescaling of the y-axis required?
        if len(Y_label) == 3:
            y_plot = [y_data*Y_label[2] for y_data in y_plot]
 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            x_plot = [x_data*X_label[2] for x_data in x_plot]
            
        returnval = ax.pcolor(x_plot, y_plot, z_plot, label=plot[3], **userargs)
        
        return ax, x_plot, returnval

    # get collection from Plot
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label, funcReturn=True, **kwargs)

# "Normal" plot
def Linear_Plot(ax, Plot_list, X_label, Y_label, **kwargs):

    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[XData, YData, "Label"],
                                 [XData2, YData2, "Label", 'linestyle=dashed'],
                                 [XData2, YData2, "Label", 'linestyle=dashed, color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """       
                
        # check dimension of X-Axis if whole trace
        x_plot = plot[0]
        
        # only one marker?
        if np.size(x_plot) > 1:
            y_plot = plot[1][0 : np.size(x_plot)]
        else:
            y_plot = plot[1]
        
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:
            
            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue
                
        # rescaling of the y-axis required?
        if len(Y_label) == 3:
            y_plot = [y_data*Y_label[2] for y_data in y_plot]
 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            x_plot = [x_data*X_label[2] for x_data in x_plot]
            
        ax.plot(x_plot, y_plot, label=plot[2], **userargs)
        
        return ax, x_plot

    # call function and return
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label, **kwargs)

# SemiLogX Plot
def SemiLogX_Plot(ax, Plot_list, X_label, Y_label, **kwargs):
    
    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[XData, YData, "Label"],
                                 [XData2, YData2, "Label", 'linestyle=dashed'],
                                 [XData2, YData2, "Label", 'linestyle=dashed, color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """       
                
        # check dimension of X-Axis if whole trace
        x_plot = plot[0]
        
        # only one marker?
        if np.size(x_plot) > 1:
            y_plot = plot[1][0 : np.size(x_plot)]
        else:
            y_plot = plot[1]
        
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:
            
            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue
                
        # rescaling of the y-axis required?
        if len(Y_label) == 3:
            y_plot = [y_data*Y_label[2] for y_data in y_plot]
 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            x_plot = [x_data*X_label[2] for x_data in x_plot]
            
        ax.semilogx(x_plot, y_plot, label=plot[2], **userargs)
        
        return ax, x_plot

    # call function and return
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label, **kwargs)

# SemiLogY Plot
def SemiLogY_Plot(ax, Plot_list, X_label, Y_label, **kwargs):
    
    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[XData, YData, "Label"],
                                 [XData2, YData2, "Label", 'linestyle=dashed'],
                                 [XData2, YData2, "Label", 'linestyle=dashed, color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """       
                
        # check dimension of X-Axis if whole trace
        x_plot = plot[0]
        
        # only one marker?
        if np.size(x_plot) > 1:
            y_plot = plot[1][0 : np.size(x_plot)]
        else:
            y_plot = plot[1]
        
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:
            
            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue
                
        # rescaling of the y-axis required?
        if len(Y_label) == 3:
            y_plot = [y_data*Y_label[2] for y_data in y_plot]
 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            x_plot = [x_data*X_label[2] for x_data in x_plot]
            
        ax.semilogy(x_plot, y_plot, label=plot[2], **userargs)
        
        return ax, x_plot

    # call function and return
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label, **kwargs)

# Boxplot
def Box_Plot (ax, Plot_list, X_label, Y_label, **kwargs):

    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[XData, YData],
                                 [XData2, YData2, 'color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """        
                
        # check dimension
        x_plot = plot[0]
        y_plot = plot[1]
        
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:

            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue
            
                
        # rescaling of the y-axis required?
        if len(Y_label) == 3:
            y_plot = [np.asarray(y_data)*Y_label[2] for y_data in y_plot]
 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            x_plot = [np.asarray(x_data)*X_label[2] for x_data in x_plot]
            
            
        # Color? -> Some Special Threatment
        if "color" in userargs:
            
            color = userargs["color"]
                        
            # patch artist
            userargs["patch_artist"] = True
            userargs["boxprops"] = dict(facecolor=color, color=color, alpha=0.75)
            userargs["capprops"] = dict(color=color)
            userargs["whiskerprops"] = dict(color=color)
            userargs["medianprops"] = dict(color='k')
            userargs["flierprops"] = dict(marker='x', markersize = 2, markeredgecolor=color)
            
            # remove "color" from list
            userargs.pop("color")
            
        # Box Widths
        if not("widths" in userargs):
            
            # calculate box width with 25% of min distance between X Points
            userargs["widths"] = np.mean( np.abs(x_plot-np.roll(x_plot,1)) )*0.25
                        
        returnval = ax.boxplot(y_plot, positions=x_plot, **userargs)
        
        # Add Legend to Axis
        handles, labels = ax.get_legend_handles_labels()
        handles.append(returnval["boxes"][0])
        labels.append(plot[2])
        ax.legend(handles, labels)
        
        return ax, x_plot, returnval

    
    # get collection from Plot
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label,
                        funcReturn=True, XAutolim=False, Legend=False, **kwargs)

# Histogram Plot
def Histogram_Plot(ax, Plot_list, X_label, Y_label, FreedmanDiacoins=True, **kwargs):

    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[XData, YData, "Label"],
                                 [XData2, YData2, "Label", 'linestyle=dashed'],
                                 [XData2, YData2, "Label", 'linestyle=dashed, color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """       
                
        # check dimension of X-Axis if whole trace
        x_plot = plot[0]
        
        # only one marker?
        if np.size(x_plot) > 1:
            y_plot = plot[1][0 : np.size(x_plot)]
        else:
            y_plot = plot[1]
        
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:
            
            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue

        # calulate bins using Freedman–Diaconis_rule
        if FreedmanDiacoins:
            q25, q75 = np.percentile(y_plot, [0.25, 0.75])
            bin_width = 2 * (q75 - q25) * len(y_plot) ** (-1/3)
            userargs["bins"] = round((max(y_plot) - min(y_plot)) / bin_width)
                
        # rescaling of the y-axis required?
        if len(Y_label) == 3:
            y_plot = [y_data*Y_label[2] for y_data in y_plot]
 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            x_plot = [x_data*X_label[2] for x_data in x_plot]
            
        ax.hist(x_plot, label=plot[2], **userargs)
        
        return ax, x_plot

    # call function and return
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label, XAutolim=False, **kwargs)

# Barh plot
def Bar_Plot(ax, Plot_list, X_label, Y_label, **kwargs):

    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[YData, Width, "Label"],
                                 [YData2, Width2, "Label", 'linestyle=dashed'],
                                 [YData3, Width3, "Label", 'linestyle=dashed, color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """       
                
        # check dimension of X-Axis if whole trace
        x_plot = plot[0]
        height_plot = plot[1]
                
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:
            
            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue
            
            # check if us bool
            if userargs[userarg] == "True":
                userargs[userarg] = True
                continue            
            if userargs[userarg] == "False":
                userargs[userarg] = False
                continue     
                 
        # rescaling of the x-axis required?
        if len(Y_label) == 3:
            height_plot = height_plot * Y_label[2]
            
        ax.bar(x_plot, height_plot, label=plot[2], **userargs)
        
        return ax, x_plot

    # call function and return
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label, **kwargs)

# Barh plot
def Barh_Plot(ax, Plot_list, X_label, Y_label, **kwargs):

    def Process_Plot(ax, plot, Y_label, X_label):
        
        """
        Prepares a X-Y linear plot
    
        paramters              description
        =====================  =============================================:
        ax                      plotting axis
        
        plot                    contains plotting array:
                                [[YData, Width, "Label"],
                                 [YData2, Width2, "Label", 'linestyle=dashed'],
                                 [YData3, Width3, "Label", 'linestyle=dashed, color=k'], ...]
                                
        Y_label                 label y axis
        
        X_label                 label x axis
        """       
                
        # check dimension of X-Axis if whole trace
        y_plot = plot[0]
        width_plot = plot[1]
                
        # emtpy argument list
        userargs = {}
                
        # insert plotting arguments
        if len(plot) >= 4:
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:
            
            # check if is int            
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                continue
                
            # check if is float
            if userargs[userarg].replace('.','',1).isdigit():
                userargs[userarg] = float(userargs[userarg])
                continue
            
            # check if us bool
            if userargs[userarg] == "True":
                userargs[userarg] = True
                continue            
            if userargs[userarg] == "False":
                userargs[userarg] = False
                continue     
                 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            width_plot = width_plot * X_label[2]
            
        ax.barh(y_plot, width_plot, label=plot[2], **userargs)
        
        return ax, y_plot

    # call function and return
    return Generic_Plot(Process_Plot, ax, Plot_list, X_label, Y_label, **kwargs)


#############################################################################
###         Generate Plot for Time Domain / Linear
#############################################################################
def Polar_Plot(ax, Plot_list, X_label, Y_label, Legend=True, LegendLoc=0,
                deg2rad = True, fontsize=14, TicksEng=True, XTicksLabel=None, 
                legendcol=1,fontsize_label=14, yaxis_pad=0, xaxis_pad=0, 
                BlackWhite=False, grid = True, minorgridalpha=0.25,
                majorgridalpha=0.5, **kwargs):
#############################################################################  
    """
    Prepares a X-Y linear plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    Plot_list               all X and Y Values also Labels (and matplotlib arguments)
    X_label                 X Axis Label and Unit (option: rescaling factor) (Engineering Package)
    Y_label                 Y Axis Label and Unit (option: rescaling factor)(Engineering Package)
    Legend                  (option) plot legend
    LegendLoc               (option) legend location
    TicksEng                (option) Enable Engineering Ticks
    XTicksLabel             (option) Label only ever nth tick
    deg2rad                 (option) Convert Angle from degree into radians
    fontsize                (option) Fontsize of the legend and ticks
    fontsize_label          (option) Fontsize of the axis labels
    legendcol               (option) Legend Columns
    yaxis_pad               (option) move label to y-axis (padding)
    xaxis_pad               (option) move label to x-axis (padding)    
    BlackWhite              (option) Use Black and White Preset
    grid                    (option) Use grid
    
    return type
       None  (writes directly into axis)
      
    Example:
        
        import basic_toolbox as basic
        
        ...
        
        # Prepare
        Xlabel = ["XLabel", 'V']
        Ylabel = ["YLabel", 'V']
        Plot = [[XData, YData, "Label"], [XData2, YData2, "Label", 'linestyle=dashed'],...]
        
        # or for loops just use
        Plot.append([XData3, YData3, "Label"])
        
        # Generate Plot
        plt.figure(figsize=(7.5,12))
        fig, ax1 = plt.subplots(subplot_kw={'projection': 'polar'})
        basic.Linear_Plot(ax1, Plot, Xlabel, Ylabel)  
        plt.show()
        
        # Two Y Axis one plot (similar to X-Axis)
        ax2 = ax1.twinx()
        basic.Linear_Plot(ax2, Plot, Xlabel, Ylabel, TwinX=ax1)         
   
    """        
#############################################################################   

    # BlackWhite Default Settings
    if BlackWhite:
        ax.set_prop_cycle(monochrome)
        
        
    for index in range(len(Plot_list)):
        
        plot = Plot_list[index]
        
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
            args = plot[3].strip().replace(" ", "")
            userargs = dict(e.split('=') for e in args.split(','))
            
        # Check if userargs have only numberic values
        for userarg in userargs:
            if userargs[userarg].isdigit():
                userargs[userarg] = int(userargs[userarg])
                
        # rescaling of the y-axis required?
        if len(Y_label) == 3:
            y_plot = [y_data*Y_label[2] for y_data in y_plot]
 
        # rescaling of the x-axis required?
        if len(X_label) == 3:
            x_plot = [x_data*X_label[2] for x_data in x_plot]
            
        # xaxis (angle) from degree to radians
        if deg2rad:
            x_plot = np.deg2rad(x_plot)
            
        ax.plot(x_plot, y_plot, label=plot[2], **userargs)
        
    # label
    #ax.set_ylabel(Y_label[0], labelpad=yaxis_pad)
    #ax.set_xlabel(X_label[0], labelpad=xaxis_pad)

    
    # ticks in engineering formatter
    if TicksEng:
        ax.yaxis.set_major_formatter(tck.EngFormatter(unit=Y_label[1]))
    

    # set font sizes (all)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
    
    # set font size label
    for item in ([ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(fontsize_label)

    # =================================== 
    # change XTick Label Position
    if XTicksLabel:
        
        # change visibility of each Nth tick
        for (index,label) in enumerate(ax.xaxis.get_ticklabels()):
            if index % XTicksLabel != 0:
                label.set_visible(False)


    # ===================================    
    # grid and legend
    else:
        
        if Legend:
            # legend
            ax.legend(framealpha=1, loc=LegendLoc, fontsize=fontsize, ncol=legendcol)
            
        if grid:
            
            ax.minorticks_on()
            ax.grid(which='major', alpha=majorgridalpha, linestyle='-',linewidth=1.2) 
            ax.grid(which='minor', alpha=minorgridalpha, linestyle=':', linewidth=1)

        
    #retrn
    return ax


#############################################################################
###         Generate Vertical Line with Label
#############################################################################
def Vline_Plot(ax, xValue, xLabel, yDistance=0.5, yPos='up', color='r',
               fontsize=6, horizontalalignment='center', **kwargs):
#############################################################################  
    """
    Generates Vertical Line in Plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    xValue                  Value on X-Axis
    xLabel                  Label for Line
    yDistance               Distance Factor (Y Tick Distance, default=0.25)
    yPos                    'up' or 'down' or 'center'
    color                   color of line and text (default=red)
    fontsize                fontsize of text (default=12)
    linestyle               linestyle of line (default='-')
    linewidth               linewidth
    horizontalalignment     Alignment of text (default='center')
    
    return type
       None  (writes directly into axis)
       
    Example:
        
        import basic_toolbox as basic

   
    """    

#############################################################################  
        
    # Add vertical line
    ax.axvline(x=xValue, color=color, **kwargs)

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
        
    if yPos == 'center': 
        ylimits = (ylimits[1] - ylimits[0])/2 + ylimits[0]
        
        

    # generate Text            
    ax.text(xValue, ylimits+ydistance, xLabel, color=color,
            fontsize=fontsize, horizontalalignment=horizontalalignment)  
    
    # jump back
    return
    
#############################################################################
###         Generate Vertical Line with Label
#############################################################################
def Hline_Plot(ax, yValue, yLabel, xDistance=0.4, yDistance=0, xPos='right',
               fontsize=6, verticalalignment='center', color='r',
               TextBG="", **kwargs):
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
        
        import basic_toolbox as basic

   
    """    

#############################################################################  
    # Add vertical line
    ax.axhline(y=yValue, color=color, **kwargs)
    yValueLabel = yValue + yDistance

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

    if xPos == 'center': 
        xlimits = (xlimits[1] - xlimits[0])/2 + xlimits[0]
        
    # generate Text            
    t = ax.text(xlimits+xdistance, yValueLabel, yLabel, color=color,
                fontsize=fontsize, verticalalignment=verticalalignment,
                horizontalalignment="center") 
    
    if TextBG:
        t.set_bbox(dict(facecolor=TextBG, alpha=0.6, linewidth=0))
    
    # jump back
    return

#############################################################################
###         Generate Fill Pattern in Plot
#############################################################################
def Fill_Plot(ax, XAxis, YAxis1, YAxis2 = None, XLimits = None, StickyLimits=True, **kwargs):
#############################################################################  
    """
    Generates Vertical Line in Plot

    paramters              description
    =====================  =============================================:
    ax                     plot axis
    XAxis                  X Axis
    YAxis1                 YAxis1
    YAxis2                 (optional) YAxis2
    XLimits                (optional) Range on the X Axis [start, end]
    StickyLimits           (optional) Stick to old YLim and XLim (first plot data)

    return type
       None  (writes directly into axis)
   
    """    

#############################################################################  

    # get old limits
    old_xlim = ax.get_xlim()
    old_ylim = ax.get_ylim()
    
    # limits for x axis
    if XLimits is not None:
        Xstart = np.argmin(np.abs(XAxis - XLimits[0]))
        Xstop = np.argmin(np.abs(XAxis - XLimits[1]))
        
        XAxis = XAxis[Xstart:Xstop]
        YAxis1 = YAxis1[Xstart:Xstop]
        
        if YAxis2 is not None:
            YAxis2 = YAxis2[Xstart:Xstop]
        
        
    if YAxis2 is None:
        YAxis2 = np.zeros(len(YAxis1))

    # generate Fill Pattern         
    ax.fill_between(XAxis, YAxis1, YAxis2,  **kwargs)

    # set old limits
    if StickyLimits:
        ax.set_xlim(old_xlim)
        ax.set_ylim(old_ylim)
        
    # jump back
    return    
    
    
#############################################################################
###         Generate Rectangle Plot
#############################################################################
def Rectangle_Plot(ax, xCenter, xSpan, yCenter, ySpan,
                   fullSpanY=False, StickyLimits=True, **kwargs):
#############################################################################  
    """
    Generates Vertical Line in Plot

    paramters              description
    =====================  =============================================:
    ax                      plot axis
    xCenter                 Rectangle x-axis center
    xSpan                   Rectangle x-axis span
    yCenter                 Rectangle y-axis center
    ySpan                   Rectangle y-axis span
    StickyLimits            (optional) Stick to old YLim and XLim (first plot data)

    return type
       None  (writes directly into axis)
       

   
    """    

#############################################################################  

    # get old limits
    old_xlim = ax.get_xlim()
    old_ylim = ax.get_ylim()
    
    # calculate startpoints
    xstart = xCenter - xSpan/2
    ystart = yCenter - ySpan/2
        
    if fullSpanY:
        ystart = min(old_ylim)
        ySpan = old_ylim[1] - old_ylim[0]

    # generate Rectangle
    rect = patches.Rectangle((xstart,ystart),xSpan, ySpan, **kwargs)          
    ax.add_patch(rect) 

    # set old limits
    if StickyLimits:
        ax.set_xlim(old_xlim)
        ax.set_ylim(old_ylim)
        
    # jump back
    return
    
#############################################################################
###         Align two Y-axis
#############################################################################
def Align_YXAxis(ax1, ax2, AxisType="Y", Method="NONE"):
#############################################################################    
    """
   Align two Axis on the same Grid

    paramters              description
    =====================  =============================================:
    ax1,ax2                 axis object
    AxisType                either "Y" or "X"
    Method                  NONE | AX1 = use AX1 Ticks | BOTH = fith both
    
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
         
    if Method == "BOTH":
        
        # try to align both axis using 
        if AxisType == "Y":
            # get maximum number of ticks
            max_ticks = max(len(ax1.get_yticks()),len(ax2.get_yticks()))
                
            # get axis distance between ticks
            ax1_dy = ax1.get_ybound()[1] - ax1.get_ybound()[0]   
            ax2_dy = ax2.get_ybound()[1] - ax2.get_ybound()[0]
            
        elif  AxisType == "X":
            # get maximum number of ticks
            max_ticks = max(len(ax1.get_xticks()),len(ax2.get_xticks()))
                
            # get axis distance between ticks
            ax1_dy = ax1.get_xbound()[1] - ax1.get_xbound()[0]   
            ax2_dy = ax2.get_xbound()[1] - ax2.get_xbound()[0]
         
        # Roundto Number
        ax1_roundto = Generate_RoundPoint(ax1_dy, max_ticks)
        ax2_roundto = Generate_RoundPoint(ax2_dy, max_ticks)
    
        
        if AxisType == "Y":     
            # get axis bounds and scale
            YBound_ax1 = [np.floor(ax1.get_ybound()[0]/ax1_roundto),
                          np.ceil(ax1.get_ybound()[1]/ax1_roundto)]
            YBound_ax2 = [np.floor(ax2.get_ybound()[0]/ax2_roundto),
                          np.ceil(ax2.get_ybound()[1]/ax2_roundto)]
    
        elif  AxisType == "X":
            # get axis bounds and scale
            YBound_ax1 = [np.floor(ax1.get_xbound()[0]/ax1_roundto),
                          np.ceil(ax1.get_xbound()[1]/ax1_roundto)]
            YBound_ax2 = [np.floor(ax2.get_xbound()[0]/ax2_roundto),
                          np.ceil(ax2.get_xbound()[1]/ax2_roundto)]  
            
        # get axis scaling and scale
        ax1_dy = YBound_ax1[1] - YBound_ax1[0]   
        ax2_dy = YBound_ax2[1] - YBound_ax2[0]
        
        # define starting points
        ax1_start = YBound_ax1[0]
        ax2_start = YBound_ax2[0]
         
        # Generate new Ticks for axis
        ax1_ticks_new =  Generate_newTicks(ax1_dy, ax1_start, max_ticks)
        ax2_ticks_new =  Generate_newTicks(ax2_dy, ax2_start, max_ticks)
        
        # reverse rescaling
        ax1_ticks_new = ax1_ticks_new*ax1_roundto
        ax2_ticks_new = ax2_ticks_new*ax2_roundto

    else:
        #skip
        return
##############################
        
    if AxisType == "Y":    
        # set ticks
        ax1.set_yticks(ax1_ticks_new)
        ax2.set_yticks(ax2_ticks_new)
            
    elif AxisType == "X":
        # set ticks
        ax1.set_xticks(ax1_ticks_new)
        ax2.set_xticks(ax2_ticks_new)       
  
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
        
        import basic_toolbox as basic
        
        # generate Dataset
        XData = ...
        YData = ...
        
        XPoint = ...
        
        
        # Find Point
        YPoint = basic.FindPoint_NextValue(XData, YData, XPoint)
           
    """
#############################################################################     
    # Search next minimum
    index = np.argmin(np.abs(XData-XPoint))

    return YData[index]

#############################################################################
###         Clock Extraction
#############################################################################
def Extract_DigitalClock(WaveX, WaveY,
                         edge_trigger = 'rising', trigger_val = 0.5,
                         threshold_high = 0.01, threshold_low = 0):
#############################################################################  

#############################################################################  
    
    # Digitialize Data
    Y_dig = np.array(Digitalize_Data(WaveY, 0, onlyDigi=True,
                                     threshold_high=threshold_high,
                                     threshold_low=threshold_low)[0])
 

    # choose edge mask for rising or falling edge clock mask
    if edge_trigger == 'rising':
        mask = (Y_dig[:-1] < trigger_val) & (Y_dig[1:] > trigger_val)
    elif edge_trigger == 'falling':
        mask = (Y_dig[:-1] > trigger_val) & (Y_dig[1:] < trigger_val) 
    else :
        print("edge trigger mask not valid!") 
            
    # Get Index and Timestamp
    mask_index = [index for index, value in enumerate(mask) if value] 
    mask_time = WaveX[mask_index]
    
    # Timestamps
    XClk = mask_time[1::]
    
    # calulate periode
    YClk = np.diff(mask_time)

    return [XClk, YClk]

#############################################################################
###         Digitalize Data
#############################################################################
def Digitalize_Data(data, clock, chipselect = [], edge_trigger = 'rising',
                    high_val = 1, low_val = 0, trigger_val = 0.5,
                    threshold_high = 2, threshold_low = 0.8, plot = False,
                    onlyDigi = False, mute= True):
#############################################################################  
    """
    Fits an higher order order polynom function against a set of data in regime 
    around a fitting point

    paramters              description
    =====================  =============================================:
    data                    data array
    clock                   clock array
    chipselect              (optional) chip select
    edge_trigger            (optional) rising or falling
    high_val                (optional) parsed high value
    low_val                 (optional) parsed low value
    trigger_val             (optional) trigger value
    threshold_high          (optional) threshold high value
    threshold_low           (optional) threshold low value
    plot                    (optional) plot data for debug
    onlyDigi                (optional) returns only data as logical data
    mute                    (optional) no print output
    
    return type
       binary bit stream
       
    Example:
        
        import basic_toolbox as basic
        
        # generate Dataset
        data = ...
        clock = ...
                
        
        # Find Point
        stream = basic.Digitalize_Data(data, clock)
           
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
                
            elif analog_value <= threshold_low:
                digital_array.append(low_val)
                
            # out of range
            else:
               conversion_errors.append(current_index)
               digital_array.append(analog_value)
               if not(mute):
                   print("Data " + str(analog_value) + " (index = " 
                         + str(current_index) + ") could not be parsed!") 
                
            # increase index
            current_index = current_index + 1
            
        return [digital_array, conversion_errors]
    
    def BinarytoByte(BinaryString):
        
        # Convert Binary String into Bytes
        Bytes = [BinaryString[i:i+4] for i in range(0, len(BinaryString), 4)]
        Bytes = [str(hex(int(i,2))).replace('0x','') for i in Bytes]
        
        # Merge every 4 Bytes
        return [ ''.join(x) for x in zip(Bytes[0::4], Bytes[1::4], Bytes[2::4], Bytes[3::4]) ]
    
#############################################################################    
       
    # Generate Digital Data
    if onlyDigi:
        return digitialize(data, high_val=high_val, 
                           low_val=low_val,
                           threshold_high=threshold_high, 
                           threshold_low=threshold_low)
    
    # Chip Select availible?
    CS = len(chipselect) == len(data)
        
        
    # Generate digital stream
    data_dig = digitialize(data, high_val=high_val, low_val=low_val,
                           threshold_high=threshold_high, threshold_low=threshold_low)
    
    clock_dig = digitialize(clock, high_val=high_val, low_val=low_val,
                            threshold_high=threshold_high, threshold_low=threshold_low)
    
    # Check for Chip Select
    if(CS):
        chipselect_dig = digitialize(chipselect, high_val=high_val, low_val=low_val,
                                     threshold_high=threshold_high, threshold_low=threshold_low)       

    # Remove Indicies
    conversion_errors = data_dig[1] + clock_dig[1]
    if(CS): conversion_errors = conversion_errors + chipselect_dig[1]
        
    # clear faulty data from streams
    data_dig = np.delete(data_dig[0], conversion_errors)
    clock_dig = np.delete(clock_dig[0], conversion_errors)
    if(CS): chipselect_dig = np.delete(chipselect_dig[0], conversion_errors)
    
    # rescale original data
    data = np.delete(data, conversion_errors)
    clock = np.delete(clock, conversion_errors)
 
    # data conversion
    data_dig = np.array(data_dig)
    clock_dig = np.array(clock_dig)
     
    # rising or falling edge clock mask
    mask_rising = (clock_dig[:-1] < trigger_val) & (clock_dig[1:] > trigger_val)
    mask_falling = (clock_dig[:-1] > trigger_val) & (clock_dig[1:] < trigger_val) 
 
    # rising or falling chip select
    if(CS):
        
        # data conversion
        chipselect_dig = np.array(chipselect_dig)
        
        CS_mask_rising = (chipselect_dig[:-1] < trigger_val) & (chipselect_dig[1:] > trigger_val)
        CS_mask_falling = (chipselect_dig[:-1] > trigger_val) & (chipselect_dig[1:] < trigger_val) 
    
        # chip select mask
        mask_cs_rising = [index for index, value in enumerate(CS_mask_rising) if value]
        mask_cs_falling = [index for index, value in enumerate(CS_mask_falling) if value]
        
        # limit to smaller length
        max_length = np.min([len(mask_cs_rising),len(mask_cs_falling)])
        mask_cs_rising = mask_cs_rising[:max_length]
        mask_cs_falling = mask_cs_falling[:max_length]
        
    # choose edge mask
    mask = []
    if edge_trigger == 'rising':
        clockmask = [index for index, value in enumerate(mask_rising) if value]
    elif edge_trigger == 'falling':
        clockmask = [index for index, value in enumerate(mask_falling) if value]
    else :
        if not(mute):
            print("edge trigger mask not valid!")         

    # Allow only data, where chip is selected
    if(CS): 
        
        binary_data = []
        binary_cs = []
        mask = []
        
        # iterate each chip select cycle
        for cs_index in range(max_length):
            
            start_index = mask_cs_rising[cs_index]
            clk_start = list(filter(lambda i: i >= start_index, clockmask))
            
            # no index found -> skip
            if (len(clk_start) == 0): 
                #print("START: {} CLK: {}".format(start_index, clk_start))
                clk_start = np.flip(list(filter(lambda i: i <= start_index, clockmask)))
                #print("SMALLER START: {} NEW CLK START {}".format(start_index, clk_start))
            
            stop_index = mask_cs_falling[cs_index]
            clk_stop = list(filter(lambda i: i >= stop_index, clockmask))  
            
            # no index found -> skip
            if (len(clk_stop) == 0):
                #print("STOP: {} CLK: {}".format(stop_index, clk_start))
                clk_stop = np.flip(list(filter(lambda i: i <= stop_index, clockmask)))
                #print("SMALLER STOP: {} NEW CLK STOP {}".format(stop_index, clk_stop))
              
            # No CLk Found
            if (len(clk_stop) == 0) and (len(clk_stop) == 0):
                #print("NO CLK FOUND IN THIS REGION")
                continue

            # get clock mask for this window
            clk_window_start = clockmask.index( clk_start[0]  )
            clk_window_stop = clockmask.index( clk_stop[0] )
            
            # clock to sample data in window
            clock_windowed = clockmask[clk_window_start : clk_window_stop]

            # append to output
            mask.append(clock_windowed )
            binary_data.append( data_dig[clock_windowed] )
            binary_cs.append( chipselect_dig[clock_windowed] ) 
            
            # plot both streams
            if plot:
                 
                # orginal data, normalized
                norm_origdata = data[start_index:stop_index] / np.max(data[start_index:stop_index])
                norm_origclock = clock[start_index:stop_index] / np.max(clock[start_index:stop_index])
                
                norm_sample = np.array(clock_windowed)- start_index
                
                plt.figure(figsize=(15,10))
                ax1 = plt.subplot(211)
                ax1.plot(data_dig[start_index:stop_index], label="Data")
                ax1.plot(norm_origdata, label="Data orig.", linewidth=0.1, color='k')
                ax1.plot(norm_sample, data_dig[clock_windowed], 'rx', label="Sampled", markersize=10)
                
                ax1.vlines(norm_sample ,0,1, color='C1', label="Clk Edge", linestyles=':')
                
                ax1.grid()
                ax1.legend()
                
                ax2 = plt.subplot(212)
                ax2.plot(clock_dig[start_index:stop_index], label="Clock")
                ax2.plot(norm_origclock, label="Clock orig" , linewidth=0.1, color='k') 
                
                ax2.vlines(norm_sample ,0,1, color='C1', label="Clk Edge", linestyles=':')
                
                ax2.grid()
                ax2.legend()
                 
                plt.show              
            
            
    else:
        
        # sample data at edge
        binary_data = data_dig[clockmask]
        mask = clockmask
        
        # plot both streams
        if plot:
            
            plt.figure(figsize=(15,10))
            ax1 = plt.subplot(211)
            ax1.plot(data_dig, label="Data")
            ax1.plot(data/np.max(data), label="Data orig", linewidth=0.1, color='k')
            ax1.plot(mask,binary_data, 'rx', label="Sampled", markersize=10)
            ax1.grid()
            ax1.legend()
            
            ax2 = plt.subplot(212)
            ax2.plot(clock_dig, label="Clock") 
            ax2.plot(clock/np.max(clock), label="Clock orig", linewidth=0.1, color='k')   
            ax2.grid()
            ax2.legend()
        
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
        
        import basic_toolbox as basic
        
        # generate Dataset
        data = ...               
        
        # Generate multiple formats 
        data = basic.CMPLX2Polar(data)
           
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
        
        import basic_toolbox as basic
        
        # generate Dataset
        dataY = ...               
        dataX = ...       
          
        # Generate multiple formats 
        dataAVG = basic.Average(dataY, dataX)
           
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
##          Frequency Domain Filtering
#############################################################################
def FrequencyFiltering(Data, pass_start, pass_end, filtertype="Bandstop", plot=False):
############################################################################# 
    """
    Average of YData, both data arrays will be splitted into eqal chunks

    paramters              description
    =====================  =============================================:
    file                    input data, which should be averaged
    
    return type
       list with data
       
    Example:
        
        import basic_toolbox as basic
        
        # generate Dataset
        file = ...  
             
        # Generate multiple formats 
        list = basic.String2List(file)
           
    """   
#############################################################################    
    # Number of Samples (Real Part)
    N_FFT =  np.int(np.size(Data)/2)
    
    # FFT of the Dataset
    Data_FFT= np.fft.fft(Data)[0:N_FFT]
    
    if pass_end == -1:
        pass_end = N_FFT
    

    # Generate Bandpassfilter
    if filtertype == "Bandpass":
            Filter_FFT = np.zeros(N_FFT)
            Filter_FFT[pass_start:pass_end] = np.ones(pass_end-pass_start)
 
    # Generate Bandstopfilter
    if filtertype == "Bandstop":
            Filter_FFT = np.ones(N_FFT)
            Filter_FFT[pass_start:pass_end] = np.zeros(pass_end-pass_start)
    
    # Plot Output
    if plot==True:
        fig, ax1 = plt.subplots(figsize=(10,5))
        ax1.plot(np.abs(Data_FFT))
        ax1.set_xlim([0, 100e3])
        ax2 = ax1.twinx()
        ax2.plot(Filter_FFT, 'k:')
                 
        plt.show()
        
    # Use Filter on Dataset 
    Data_Filtered = Data_FFT * Filter_FFT
    
    # Generate Real Spektrum
    Data_Filtered = np.append(Data_Filtered, np.flip(Data_Filtered))
    
    # Return IFFT
    return np.fft.ifft(Data_Filtered)
    
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
        
        import basic_toolbox as basic
        
        # generate Dataset
        file = ...  
             
        # Generate multiple formats 
        list = basic.String2List(file)
           
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
  
#############################################################################
##          XYZ-Plot
#############################################################################
def XYZ_Plot(raw_data, xname="Time[s]", yname="Freq[Hz]", zname="Mag[dBm]", nan=float("NaN")):
############################################################################# 
    """
    Average of YData, both data arrays will be splitted into eqal chunks

    paramters              description
    =====================  =============================================:
    raw_data                input dict with 3 levels
    xname                   name of x dict
    yname                   name of y dict
    zname                   name of z dict
    
    return type
       x,y and grid matrix
                  
    """   
#############################################################################    
    data_xyz = {}
    grid_y = []
    grid_x = raw_data[xname]
    
    # Generate Dictionary with Values
    for xindex in range(len(raw_data[xname])):
        
        # List with equal time
        x = raw_data[xname][xindex]
        same_x = {}
        
        for yindex in range(len(raw_data[yname][xindex])):
                
            # extract data
            y = raw_data[yname][xindex][yindex]
            z = raw_data[zname][xindex][yindex]
            
            # generate y axis
            if y not in grid_y:
                grid_y.append(y)
    
            # insert into array
            same_x[y] = z
            
        data_xyz[x] = same_x
    
    # sort list by value
    grid_y.sort()
    
    # empty data matrix    
    grid_xyz = np.zeros((len(grid_y),len(grid_x)))

############################################################################# 
    
    # fill matrix with points
    for index_x in range(len(grid_x)):
        for index_y in range(len(grid_y)):
            
            # can be parsed?
            try:
                z_value = data_xyz[grid_x[index_x]][grid_y[index_y]]
            except:
                z_value = nan
                
            grid_xyz[index_y][index_x] = z_value
        

    return [grid_x, grid_y, grid_xyz]   


#############################################################################
##          Spectrum_Minimizer and Average
#############################################################################
def Spectrum_Minimizer(Freq_Matrix, Mag_Matrix, nanvalue=-100, minmax=True,
                       Freqsize=30, Freq_min=0, Freq_max=1, plotall=False,
                       plotnonlock=False, diff_max_lock=40, plotname=""):
############################################################################# 
    """
    Minimizes Spectrum fror Small Short Plots

    paramters              description
    =====================  =============================================:
    Freq_Matrix             Matrix with Frequency Values in Hz
    Mag_Matrix              Matrix with Magnitue Values in dBm
    nanvalue                value of NaN
    minmax                  (optional)
    
    return type
       Freq_New_Array, Maximum_Spectrum, Mean_Spectrum as Dict
                  
    """   
#############################################################################      

    # check size of both matricies
    if not(Freq_Matrix.shape == Mag_Matrix.shape):
        print("ERROR: Size of input matrices not equal!")

    # find min and max value
    if minmax == False:
        Freq_min = np.min(Freq_Matrix)
        Freq_max = np.max(Freq_Matrix)    
        
    # Generate new Freq Axis, equal spacing    
    Freq_New_Array = np.linspace(Freq_min,Freq_max, Freqsize)

    # length of samples
    samples = Freq_Matrix.shape[0]
        
    # empty resulting matrix with NaN
    Matrix = np.zeros([samples, Freqsize])
    Matrix[:] = nanvalue

    # iterate all entries in freq matrix
    for index_x in range(Freq_Matrix.shape[0]):
        for index_y in range(Freq_Matrix.shape[1]): 
            
            # extract values from given matrix
            freq = Freq_Matrix[index_x][index_y]
            mag = Mag_Matrix[index_x][index_y]
            
            # matching freq to new matrix
            match_index = np.argmin(np.abs(Freq_New_Array-freq))
            
            if not(Matrix[index_x][match_index] == nanvalue):
                old_value = Matrix[index_x][match_index]
                Matrix[index_x][match_index] = np.mean([mag,old_value])
            else:
                Matrix[index_x][match_index] = mag
    

#############################################################################  
      
    # Spectrum
    Maximum_Spectrum = []
    Minimum_Spectrum = []
    Mean_Spectrum = []
    
    Maximum_onlyValues = []
    
    # generate Max and Meal Value
    for column in Matrix.T:
        
        # find all nanvalues from 
        nanvalue_index = []
        for index, item in enumerate(column):
            if item == nanvalue:
                nanvalue_index.append(index)
             
        # remove all nanvalues
        column = np.delete(column, nanvalue_index)
        
        # generate ouput
        if column.any():
            Maximum_Spectrum.append(np.max(column))
            Minimum_Spectrum.append(np.min(column))
            Maximum_onlyValues.append(np.max(column))
            Mean_Spectrum.append(np.mean(column))
            
        else:
            Maximum_Spectrum.append(float("NaN"))
            Minimum_Spectrum.append(float("NaN"))
            Mean_Spectrum.append(float("NaN"))           
     

#############################################################################         
    # generate peak and mean
    mean_max = np.mean(Maximum_onlyValues)
    peak_max = np.max(Maximum_onlyValues)
    
    # difference
    diff_max = peak_max-mean_max
    
    text = "diff_max: " + str(diff_max)
    
    if diff_max < diff_max_lock:
        lock = False
    else:
        lock = True
    
#############################################################################  
    
    if plotall or (plotnonlock and not(lock)):
        
        # plot the data
        plt.figure(figsize=(6,5))
        plt.plot(Freq_New_Array, Maximum_Spectrum)
        plt.plot(Freq_New_Array, Mean_Spectrum)
        plt.plot(Freq_New_Array, Minimum_Spectrum) 
        plt.title(plotname)
        
        x = np.mean(Freq_New_Array)
        y = mean_max
        plt.text(x, y, text, horizontalalignment='center',
                 verticalalignment='center')
        
        plt.grid()
        plt.show()
         
#############################################################################  
	
	# Return type
    return_dict = {}
    return_dict["Freq[Hz]"] = Freq_New_Array
    return_dict["Max_Spectrum[dBm]"] = Maximum_Spectrum
    return_dict["Mean_Spectrum[dBm]"] = Mean_Spectrum
    return_dict["Lock_Estimation"] = lock
    
    return return_dict


#############################################################################
##          Moving Filter
#############################################################################
def MovingFilter(Xdata, Ydata, FilterType="MovingAvg", **kwargs):
############################################################################# 
    """
    Simple Moving Filter

    paramters              description
    =====================  =============================================:
    Ydata                   Filter Y Data
    Xdata                   Correspnding X Data
    FilterType              Select Filter Type
    *kwargs                 Additional Filter Parameters                  
    
    return type
       YData, XData as Dict
                 
    """  
    # https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    # moving average
    def MovingAverage(Xdata, Ydata, N=3):
        
        cumsum, moving_aves = [0], []
        
        for i, x in enumerate(Ydata, 1):
            cumsum.append(cumsum[i-1] + x)
            if i >= N:
                moving_ave = (cumsum[i] - cumsum[i-N])/N
                moving_aves.append(moving_ave)
                
        return_dict = {}
        return_dict["YData"] = np.array( moving_aves )        
        return_dict["XData"] = np.array( Xdata[int((N-1)/2):-1*int((N-1)/2)] )
        
        return return_dict
 
    # ==============================================================
    # https://stackoverflow.com/questions/34235530/how-to-get-high-and-low-envelope-of-a-signal
    # envelope
    def Envelope(Xdata, Ydata, dmin=1, dmax=1, split=False):
        """
        Input :
        s:          1d-array, data signal from which to extract high and low envelopes
        dmin, dmax: int, optional, size of chunks, use this if the size of the 
                    input signal is too big
        split:      bool, optional, if True, split the signal in half along its mean
                    might help to generate the envelope in some cases
        
        """
    
        # locals min and maxima     
        lmin = (np.diff(np.sign(np.diff(Ydata))) > 0).nonzero()[0] + 1 
        lmax = (np.diff(np.sign(np.diff(Ydata))) < 0).nonzero()[0] + 1 
        
        if split:
            # s_mid is zero if s centered around x-axis or more generally mean of signal
            s_mid = np.mean( Ydata )
            
            # pre-sorting of locals min based on relative position with respect to s_mid 
            lmin = lmin[ Ydata[lmin]<s_mid]
            
            # pre-sorting of local max based on relative position with respect to s_mid 
            lmax = lmax[ Ydata[lmax]>s_mid]
    
        # global min of dmin-chunks of locals min 
        lmin = lmin[[i+np.argmin(Ydata[lmin[i:i+dmin]]) for i in range(0,len(lmin),dmin)]]
        
        # global max of dmax-chunks of locals max 
        lmax = lmax[[i+np.argmax(Ydata[lmax[i:i+dmax]]) for i in range(0,len(lmax),dmax)]]

        return_dict = {}
        
        return_dict["YData_min"] = Ydata[lmin] 
        return_dict["XData_min"] = Xdata[lmin]
        
        return_dict["YData_max"] = Ydata[lmax] 
        return_dict["XData_max"] = Xdata[lmax]

        return return_dict
  
    
    # ===============================================================
	# Select Filter

    if FilterType == "MovingAvg":
        return MovingAverage(Xdata, Ydata, **kwargs)
        
    if FilterType == "Envelope":
        return Envelope(Xdata, Ydata, **kwargs)
          


#############################################################################
##          Limit Dataset using Center
#############################################################################
def LimitDataset(Xdata, Ydata, Xcenter, Xrange):
############################################################################# 
    """
    Simple Moving Filter

    paramters              description
    =====================  =============================================:
    Ydata                   Filter Y Data
    Xdata                   Correspnding X Data
    N                       (optional) Moving Filter Size
    minmax                  
    
    return type
       YData, XData as Dict
                 
    """  
   
    # Find Center Index
    XCenter_index = np.argmin( np.abs(Xdata - Xcenter) )
     
    # return List
    return_dict = [] 
    return_dict.append( Xdata[ XCenter_index-Xrange : XCenter_index+Xrange ] )
    return_dict.append( Ydata[ XCenter_index-Xrange : XCenter_index+Xrange ] )    
    return return_dict 

#############################################################################
###        Print Progressbar
#############################################################################
def printProgressBar (iteration, total, prefix = "Progress", suffix = "Complete",
                      decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    Based on: https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
  
 #############################################################################
 ###        Engineering Notation
 #############################################################################       
def EngNot( x, sig_figs=3, si=True, exp=True):
    """
    Returns float/int value <x> formatted in a simplified engineering format -
    using an exponent that is a multiple of 3.

    sig_figs: number of significant figures

    si: if true, use SI suffix for exponent, e.g. k instead of e3, n instead of
    e-9 etc.
    Based on https://stackoverflow.com/questions/17973278/python-decimal-engineering-notation-for-mili-10e-3-and-micro-10e-6
    """
    
    if np.isnan(x):
        return ""
    
    
    x = float(x)
    sign = ''
    if x < 0:
        x = -x
        sign = '-'
    if x == 0:
        exp = 0
        exp3 = 0
        x3 = 0
    else:
        exp = int(np.floor(np.log10( x )))
        exp3 = exp - ( exp % 3)
        x3 = x / ( 10 ** exp3)
        x3 = round( x3, -int( np.floor(np.log10( x3 )) - (sig_figs-1)) )
        if x3 == int(x3): # prevent from displaying .0
            x3 = int(x3)

    if si and exp3 >= -24 and exp3 <= 24 and exp3 != 0:
        exp3_text = 'yzafpnum kMGTPEZY'[ exp3 // 3 + 8]
    elif exp3 == 0:
        exp3_text = ''
    else:
        exp3_text = 'e%s' % exp3

    return ( '%s%s%s') % ( sign, x3, exp3_text)
  
 #############################################################################
 ###        Reverse Engineering Notation
 #############################################################################       
def RevEngNot( x ):
    """
    Returns float/int value <x> from SI formatted input 
    """
    
    pos_postfixes = ['k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
    neg_postfixes = ['m', 'u', 'n', 'p', 'f', 'a', 'z', 'y']

    num_postfix = x[-1]
    
    if num_postfix in pos_postfixes:
        num = float(x[:-1])
        num*=10**((pos_postfixes.index(num_postfix)+1)*3)
    elif num_postfix in neg_postfixes:
        num = float(x[:-1])
        num*=10**(-(neg_postfixes.index(num_postfix)+1)*3)
    else:
        num = float(x)
 
    if np.isnan(num):
        return ""
             
    return num
