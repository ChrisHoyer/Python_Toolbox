#############################################################################
#   Verschiedene Skripte um ADS Simulationen und Mappings
#   - ImportLogfile
#   - ImportMapping    
#   - ImportDifferentialMapping
#   - Mapping2DiffMatrix
#   - FindWorstCells
#
#   Autor: C. Hoyer (choyer.ch@gmail.com)
#   Stand: 19.03.19
#############################################################################

import skrf as rf
import numpy as np
import re
import csv
from os import listdir

#############################################################################
###         Import Mapping Logile
#############################################################################
def ImportLogfile(file, subrow='DCBIAS',  **kwarg):
#############################################################################    
    """
   Import Logfile from Mapping

    paramters              description
    =====================  =============================================:
    file                    filepath to .log
    subrows                 interrested row (eg. DCBIA or SPARAMETER)
    
    return type
       dict for each cell
    
    Example:
        
        import Mapping_Toolbox as mapping
        
        log_file = '../../Schaltungen/ABB004GH020_CHO/Mapping_CurrentSource/biaslist/C1807-3_ABB004GH_20_CHO_current.log'
 
        Logfile = mapping.ImportLogfile(log_file)
        
        
    """
#############################################################################        
    # Import TXT Logfile
    file = open(file, 'r')
        
    # generate new Block
    Block = {}
                 
    # Import Block
    for row in file:
         
        # compress data
        compressrow = row.replace(":", "")
        compressrow = compressrow.replace(" ", "") 
        
        # Split String with Tab
        subrows = re.split(r'\t+', compressrow)
        
        # Skip if there is only one Tab
        if len(subrows) == 1:
            continue
    
        # seperate for Cellname
        cellname = subrows[1]   
         
        # check dictionary
        if not(cellname in Block.keys()):
            Block[cellname] = {} 
            
        rowname = subrows[2]
        
        # check Data
        if rowname == subrow:
             
            # Get value
            rowdata =  subrows[3] 
            rowdata = rowdata.split(',')
            
            # Get all Data
            BlockData = Block[cellname]
            
            # iterate all data
            for datapair in rowdata:
                
                # split data
                pair = datapair.split('=')
                pairkey = pair[0]
                pairdata = pair[1]            
                
                # check dictionary
                if not(pairkey in BlockData.keys()):
                        BlockData[pairkey] = [] 
                
                # convert Data
                numeric = '0123456789-.'
                for i,c in enumerate(pairdata):
                    if c not in numeric:
                        break
                    
                pairvalue = float(pairdata[:i])
                pairunit = pairdata[i:].lstrip()  
                
                # Seperate V, A, W
                pairunit = pairunit.replace(" ", "")
                pairunit = pairunit.replace("V", "")
                pairunit = pairunit.replace("A", "")
                pairunit = pairunit.replace("W", "")
                
                # Check SI-prefix
                if pairunit == "m":
                    pairvalue = pairvalue*1e-3
                
                #print(str(pairvalue))
                
                # Add Data to BlockData
                BlockData[pairkey].append(pairvalue)
               
            #Add Blockdata to Block
            Block[cellname] = BlockData
    
    return Block
    
#############################################################################
###         Import  Mapping from Folder
#############################################################################
def ImportMapping(folder, ZParam = False, **kwargs):
#############################################################################    
    """
   Import Mapping from  folder

    paramters              description
    =====================  =============================================:
    folder                 filepath to mapping folder
    ZParam                 convert to Z paramter
    
    return type
       dict for each cell
  
    """
#############################################################################  
    # Export File
    Mapping = {}
    
    # Look for s2p files in folder
    s2p_files = []
    for file in listdir(folder):
        if file.endswith(".s2p"):
            s2p_files.append(file)
            
    # Get CircuitName
    CircuitName = s2p_files[-1].replace('_c', ' ')
    CircuitName = CircuitName.split()[0]
    
    # Info
    print("============")
    print("Import Ciruit: " + CircuitName + " with " + str(len(s2p_files)) + " Cells")  
    
    # Get Cell Names
    for Cell in s2p_files:
        Cellname = Cell.replace(CircuitName, '')
        Cellname = Cellname.replace('.s2p', '')
        Cellname = Cellname.replace('_', '')
        
        # generate sweeps
        s2p_cell_file = folder + Cell
          
        # Import S2P matrix
        s2p_matrix = rf.Network(s2p_cell_file)
  
        # Empty Matrix
        Block_Sdd = {}
        Block_Zdd = {} 
        
        # Differential Matrix
        Block_Sdd['S11'] = s2p_matrix.s11.s.reshape(-2)  
        Block_Sdd['S12'] = s2p_matrix.s12.s.reshape(-2)   
        Block_Sdd['S21'] = s2p_matrix.s21.s.reshape(-2)  
        Block_Sdd['S22'] = s2p_matrix.s22.s.reshape(-2)  
        
        # Import Z-Matrix
        if ZParam:
            Block_Zdd['Z11'] = s2p_matrix.z[:,0,0]
            Block_Zdd['Z12'] = s2p_matrix.z[:,0,1]
            Block_Zdd['Z21'] = s2p_matrix.z[:,1,0]
            Block_Zdd['Z22'] = s2p_matrix.z[:,1,1]

        Data = {}
            
        # Generate DB and Phase
        for Cleantag in Block_Sdd:
            Data[Cleantag] = Block_Sdd[Cleantag]
            Data[Cleantag + "_dB"] = 20*np.log10(np.absolute(Block_Sdd[Cleantag]))
            Data[Cleantag + "_Angle"] = np.angle(Block_Sdd[Cleantag]) 
            
            # extract phase in degree and unwrapped
            phase = np.angle(Block_Sdd[Cleantag])
            phase = np.unwrap(phase)
            Data[Cleantag + "_Phase"] = phase*360/(2*np.pi)
            
        # Process Z-Matrix
        if ZParam:
             for Cleantag in Block_Zdd:
                 Data[Cleantag] = Block_Zdd[Cleantag]
                 Data[Cleantag + "_mag"] = np.absolute(Block_Zdd[Cleantag])               
                
        # frequency
        freq = s2p_matrix.frequency.f
        Data['freq'] = freq
        
        # import in fstructure
        Mapping[Cellname] = Data
    
    print("Import of Mapping complete") 
    
    return Mapping
    
#############################################################################
###         Import Differential Mapping from Folder
#############################################################################
def ImportDifferentialMapping(folder, Trace_folder = 'TRACE/', ZParam=False, Z0=50, **kwargs):
#############################################################################    
    """
   Import Mapping from  folder

    paramters              description
    =====================  =============================================:
    folder                 filepath to mapping folder
    Trace_folder           subfolder for more parameters
    ZParam                 (option) calculate Z Parameters
    Z0                     (option) charateristic impedance
    
    return type
       dict for each cell
  
    """
############################################################################# 
    # Export File
    Mapping = {}
    
    # Look for s2p files in folder
    s2p_files = []
    for file in listdir(folder):
        if file.endswith(".s2p"):
            s2p_files.append(file)
            
    # Get CircuitName
    CircuitName = s2p_files[-1].replace('_c', ' ')
    CircuitName = CircuitName.split()[0]
    
    # Info
    print("============")
    print("Import Ciruit: " + CircuitName + " with " + str(len(s2p_files)) + " Cells")  
    
    # Get Cell Names
    for Cell in s2p_files:
        Cellname = Cell.replace(CircuitName, '')
        Cellname = Cellname.replace('.s2p', '')
        Cellname = Cellname.replace('_', '')
        
        # generate sweeps
        s2p_cell_file = folder + CircuitName + '_' + Cellname + '.s2p'
        swp_cell_file = folder + Trace_folder + CircuitName + '_' + Cellname + '.swp'
          
        # Calulate Mapping
        Data = Mapping2DiffMatrix(s2p_cell_file, swp_cell_file, ZParam=ZParam, Z0=Z0)
        
        # import in fstructure
        Mapping[Cellname] = Data
    
    
    print("Import of differential Mapping complete") 
    
    return Mapping
    
#############################################################################
###         Extract Differential Matrix from mapping (swp + s2*)
#############################################################################
def Mapping2DiffMatrix(diff_file, sweep_file, freqscaling=1.0e9, ZParam=False, Z0=50, **kwargs):
#############################################################################    
    """
   Import Mapping from  folder

    paramters              description
    =====================  =============================================:
    diff_file               filepath to s2p file
    sweep_file              filepath to swp file
    freqscaling            (option) freqscaling
    ZParam                 (option) calculate Z Parameters
    Z0                     (option) charateristic impedance    
    return type
       dict for each cell with 4x4 Sparam Matrix
  
    """
#############################################################################        
    # Import S2P matrix
    S2P_diff = rf.Network(diff_file)
    
#############################################################################    
    # Import SWP Matrix
    csv_reader = csv.reader(open(sweep_file),  dialect="excel-tab")
    
    # generate new Block
    line_count = 0
    Block = {}
             
    # Import Block
    for row in csv_reader:
    
        #print(f'{row}') 
        
        # Skip Global Header
        if  any('!' in s for s in row):
            continue
        
        # Find beginning of the File
        if line_count == 0:
            header = [r.replace(' ', '') for r in row]
            
            # Generate Columnheader
            for h in header:
                Block[h] = [] 
                
        # skip second line with datatype
        if line_count == 1:
            line_count += 1
            continue
        
        # get Data
        if line_count > 1:
            
            # extract data
            for h, v in zip(header, row):
                # Convert Data
                try:
                    v = float(v)
                except:
                    v = 0
                    
                # Insert Data
                Block[h].append(v)
                
        # next line    
        line_count += 1
#############################################################################
        
    Data = {}
    
    # Frequency row
    Data['freq'] = np.asarray(Block['Freq'])*float(freqscaling)
    
        
    # get separate Headers
    Re_Header = [s for s in header if "Re" in s] 
    
    # Get Clean Headername
    for header in Re_Header:
        Cleantag = header.replace('Re', '')
        
        # generate complex Value
        Realvalue = np.asarray(Block['Re'+Cleantag])
        Imagvalue = np.asarray(Block['Im'+Cleantag]) 
        Complex =  Realvalue + Imagvalue*1j
        
        # generate dB and Phase
        Data[Cleantag] = Complex
        Data[Cleantag + "_dB"] = 20*np.log10(np.absolute(Complex))
        Data[Cleantag + "_Phase"] = np.angle(Complex)  
            
#############################################################################   
    
    # Differential Matrix
    Block_Sdd = {}
    Block_Sdd['Sdd11'] = S2P_diff.s11.s.reshape(-2)  
    Block_Sdd['Sdd12'] = S2P_diff.s12.s.reshape(-2)   
    Block_Sdd['Sdd21'] = S2P_diff.s21.s.reshape(-2)  
    Block_Sdd['Sdd22'] = S2P_diff.s22.s.reshape(-2)   
 
    # Generate DB and Phase
    for Cleantag in Block_Sdd:
        Data[Cleantag] = Block_Sdd[Cleantag]
        Data[Cleantag + "_dB"] = 20*np.log10(np.absolute(Block_Sdd[Cleantag]))
        Data[Cleantag + "_Phase"] = np.angle(Block_Sdd[Cleantag], deg=True)  

#############################################################################   
    if ZParam:
        
        # impedance
        Data['Z0'] = S2P_diff.z0 

        # index of iteration
        index = 0
        
        # generate empty Export Matrix
        ZMatrix = {}
        index_column = [1,2,3,4]
        index_row =  [1,2,3,4]
        for row in index_row:
            for col in index_column:
                ZMatrix['Z' + str(col) + str(row)] = []
                    
        # iterate each matrix
        for Z0 in Data['Z0']:
            
            # generate some useful arrays
            M_ident = np.identity(4)
            M_sqrtZ0 = np.identity(4)*np.sqrt(Z0[0])
            
            # get S Param Data (Kann man bestimmt eleganter lösen..)
            S1x = [Data['Sdd11'][index], Data['Sdd12'][index], Data['Sdc11'][index], Data['Sdc12'][index]]
            S2x = [Data['Sdd21'][index], Data['Sdd22'][index], Data['Sdc21'][index], Data['Sdc22'][index]]
            S3x = [Data['Scd11'][index], Data['Scd12'][index], Data['Scc11'][index], Data['Sdc12'][index]]
            S4x = [Data['Scd21'][index], Data['Scd22'][index], Data['Scc21'][index], Data['Sdc22'][index]]
            
            # Generate SMatrix
            Sparam_Matrix = [S1x,S2x,S3x,S4x]
            
            # Generate ZMatrix
            Z_Invers = np.linalg.inv(np.matrix(M_ident-Sparam_Matrix))
            Z_Normal = np.matrix((M_ident+Sparam_Matrix))
            
            Zparam_Matrix = np.dot(Z_Invers,M_sqrtZ0)
            Zparam_Matrix = np.dot(Z_Normal,Zparam_Matrix)
            Zparam_Matrix = np.dot(M_sqrtZ0,Zparam_Matrix)
            Zparam_Matrix = np.asarray(Zparam_Matrix)
            
            # fill Std and MixedMode Matrix   
            for row in index_row:
                for col in index_column:
                    ZMatrix['Z' + str(col) + str(row)].append(Zparam_Matrix[col-1][row-1])
                    
            #index++
            index = index + 1
            
        # Combine to return typ
        for Cleantag in ZMatrix:
            Data[Cleantag] = ZMatrix[Cleantag]
            Data[Cleantag +'_mag'] = np.abs(ZMatrix[Cleantag])
    
    
#############################################################################    
    return Data    

#############################################################################
###         Find worst Cell numbers (Mixed Mode!)
#############################################################################   
def FindWorstCells(Mapping, Mapping_Cell, Simulation, Simulation_Cell, mean_Tolerance,
                   differential = False, debug=False, tendencycheck = -1.0, **kwargs):
#############################################################################    
    """
    Search for cells that are not in tolerance band around simulation

    paramters              description
    =====================  =============================================:
    Mapping                 Mapping Matrix
    Mapping_Cell            Compare this parameter from mapping
    Simulation              Simulation Matrix
    Simulation_Cell         Compare this parameter from Simulation
    mean_Tolerance          Allowed Tolerance
    differential            (option) Differencial Matrix?
    debug                   (option) Show all results, debuginfo
    tendencycheck           (option) allow better results
    
    return type
       dict for each cell with 4x4 Sparam Matrix
       
    Example:
        
        import Mapping_Toolbox as mapping
        
        # Import Mapping
        folder = '../../Schaltungen/ABB004GH020_CHO/Mapping/'
        Mapping = mapping.ImportMapping(folder, Z2P=True)
        
        # import Simulation
        Sim_file = '../../Schaltungen/ABB004GH020_CHO/ABB004GH020_CHO_Stability.s2p'
        Simulation = ads.ImportS2P(Sim_file, ZParam=True)
        
        Skipped_Cells = mapping.FindWorstCells(Mapping, 'S21_dB',
                                       Simulation, 'S21_dB', 4.6,
                                       debug=False, tendencycheck = 1.0)
        
          
    """
#############################################################################         
    FoundCells = []
    
    # iterate Cells
    for Cellname in Mapping:
            
        # get Cellcontent
        Cell = Mapping[Cellname]

        # Errorlist
        Error = []
        max_Error = 0
        max_Error_freq = 0
        skip = 0
        
        # Iterate all frequency components
        for freq in Cell['freq']:
                        
            # get reference frequency index
            nearest_freq = np.abs(np.asarray(Simulation['freq'])-freq)
            Sim_index = np.argmin(nearest_freq)
                        
            # get measurement frequency index
            Ref_index = np.where(Cell['freq'] == freq)
                        
            # Signals
            if differential:
                Sim_freq = Simulation['MixedMode_dB'][Simulation_Cell][Sim_index]
            else:
                Sim_freq = Simulation[Simulation_Cell][Sim_index]
                
            Ref_freq = Cell[Mapping_Cell][Ref_index]
                    
            # Search Cells out of Tolerance
            Error.append(Ref_freq-Sim_freq)
            
            # get max Error
            if np.abs(max_Error) < np.abs(Error[-1]):
                max_Error = Error[-1]
                
                   
        # Calculate mean error of cell
        mean_Error = np.mean(np.abs(Error))
        

        # check mean Error
        if (mean_Error > mean_Tolerance):
            
            if (np.sign(np.mean(Error)) == tendencycheck):
                # Add to List
                FoundCells.append(Cellname)
                skip = 1
                
        # check if lowest Error out of tolerance
        if (np.abs(Error[0])  > mean_Tolerance ):
            
            if (np.sign(Error[0]) == tendencycheck):
                # Add to List
                FoundCells.append(Cellname)
                skip = 1
                                         
                
        if debug:
            
            text = "Cell: " + Cellname + " mean Error: " + str(mean_Error) 
            text = text + " max. Error: " + str(max_Error[0]) 
            text = text + " max. Error Start: " + str(Error[0][0]) 
            text = text + " Tendency: " + str(np.sign(np.mean(Error))) 
            text = text + " Skip: " + str(skip)
            
            W  = '\033[0m'  # white (normal)
            R  = '\033[31m' # red

            # mark skipped
            if skip:
                print(R + text + W)
            else:
                print(text)
                    

    # Plot info
    print('Yield: ' + str(len(Mapping)-len(FoundCells)) + ' out of ' + str(len(Mapping)))
    if FoundCells:
        print('Skipped Cells: ' + ', '.join(map(str, FoundCells)))
               
    return FoundCells 
                