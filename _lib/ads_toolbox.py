#############################################################################
#   Verschiedene Skripte um ADS Simulationen auszuwerten
#
#   - ImportData: Imports an ASCII file from ads
#   - Calc_HarmonicBalance_Single: Calculates IP3, Psat, 1dB Comp based on HB Sim
#   - Calc_HarmonicBalance: Uses Calc_HarmonicBalance_Single for multiple frequencies
#   - Calc_StabGain: Berechnet MSG/MAG und K-Fct auf Basis von S2P
#   - Calculate_StabCircle: Berechnet Stabilitätskreise für eine Freq
#   - MixedModeSparam: Extrahiert Mixed Mode Paramters aus Simulation
#   - ImportS2P: Importiert S2P File und berechnet Z und H Parameter
#
#   Autor: C. Hoyer (info@chrishoyer.de)
#   Stand: 26.03.19
#############################################################################

import csv
import numpy as np
import skrf as rf
import matplotlib.pyplot as plt

#############################################################################
###         Import ASCII Data from ADS File
#############################################################################

def ImportData(file, Rename=False, Rename_Table = [], Merge=False, MergeColumn='', complexdata=False, **kwargs):
#############################################################################    
    """
	OBSOLTE! PLEASE LOOK AT BASIC_TOOLBOX
   Import ADS ASCII File with Data (e.g. HB Simulation)
    
    file format should look like:
        
        F0	Pin	  dBm(HB.Vout[::,1])[2, ::]
        3	-1	  1.21554116043965270E0
        3	-9	  2.20127736891134520E0
        3	-8	  3.18563038209765810E0

    parameters              description
    =====================  =============================================:
    file                    path to file
    Rename                  boolean to rename columns  
    Rename_Table            rename mapping (see example)
    Merge                   boolean to merge equal columns
    MergeColumn             name if that column
    complexdata             boolean for complex format 3.34 / -5.11  (mag + angle)
    
    return type
       dictionary with dictionary inside then listed with data 
    
    Example:
        
        import ADS_Toolbox as ads
        
        Rename_Table= {'dBm(HB.Vout[::,1])' : "Pout_1"}
        Rename_Table['dBm(HB.Vout[::,3])'] = "Pout_3"

        # Import path
        CS_2x60 = '../Kennlinien/HB/2x60_12V_1V25.txt'

        # Import Simulation
        Data = ads.ImportData(CS_2x60, Rename=True, Rename_Table=Rename_Table, Merge=True, MergeColumn='F0')
    """
#############################################################################    
    # open file
    csv_reader = csv.reader(open(file),  dialect="excel-tab")
    
    # generate new Block
    line_count = 0
    block_count = 0
    Block = {}
    Data = []
         
     # Import Block
    for row in csv_reader:
        
        # New Block if empty Line
        if not row:
           
            # add if it is not empty
            if Block:
                Data.append(Block)
            Block = {}
            block_count += 1
            line_count = 0
            continue
        
        # New Header and Content
        if line_count == 0:
            header = row
            
            # Generate Columns
            for h in header:
                Block[h] = []
            
            line_count += 1
            
        else:
            for h, v in zip(header, row):
                try:
                    
                    # check if data is complex
                    if complexdata:
                        
                        # Split at seperator
                        v = v.split('/')
                        
                        if len(v) == 2:
                            v_mag = float(v[0])
                            v_angl = float(v[1])
                        
                            # only manitude?
                            v = v_mag * np.exp(1j*np.deg2rad(v_angl))
                            
                        else:
                            v = float(v[0])
                                      
                    else:
                        # only maginute
                        v = float(v)
                    
                except:
                    
                    v = 0
                    
                Block[h].append(v)
            line_count += 1
                
    # last Block
    if Block:
        Data.append(Block)
        
#############################################################################
# Rename Datadicts
    
    if Rename:
        for datadict in Data:

            # List all Keys in this Dataset
            Dataset_Keys = list(datadict.keys())
                
            for Rename_Item in Rename_Table:

                # list of availible keys
                Dataset_newKey = []
                
                # look for Key in Dataset
                for Dataset_Key in Dataset_Keys:
                    
                    # rename found
                    if Rename_Item in Dataset_Key:
                        
                        # get new Key
                        Dataset_newKey = Rename_Table.get(Rename_Item,"");
                    
                        # item found, rename
                        if Dataset_newKey:
                            datadict[Dataset_newKey] = datadict[Dataset_Key]
                            del datadict[Dataset_Key]   
                            
 #############################################################################               
# Merge if two Datasets with same Parameters
    
    if Merge:
        
        # Looking for same Column
        MergeColum_List = []
        
        for datadict in Data:
            
            try:
                MergeColum_List.append(datadict[MergeColumn][0])
            except:
                print('Error in Merging Data, index:' + str(Data.index(datadict)))
            
        MergeColum_same = { value : [ i for i, v in enumerate(MergeColum_List) if v == value ] for value in set(MergeColum_List) }
        
        for same_Column in MergeColum_same:
                    
            index = MergeColum_same[same_Column]
           
            # Check index length
            keys_dictionary = list()
            for sub_index in range(0,(len(index))):
                    
                # first row, skip
                if (sub_index == 0):
                    continue             
                    
                # extract all new keys
                new_keys = list(Data[index[sub_index]].keys())
                    
                # Find availiable keys
                for key in set(keys_dictionary) & set(new_keys):
                    new_keys.remove(key)
                         
                # Add Content and new key to index 0
                for key in new_keys:        
                    Data[index[0]][key] = Data[index[sub_index]][key]
                    del Data[index[sub_index]][key]   
                    
    # return
    return Data

#############################################################################
###         Extract Harmonic Balance and Plot for one Freq
#############################################################################

def Calc_HarmonicBalance_Single(Data, Freq, lin_point_shift=0, lin_point_tolerance = 0.25, forcezero=False,  **kwargs):
#############################################################################   
    """
    Calculate Harmonic Balance Values for a single frequency
    
    Data format should be from the function: ImportData()

    parameters              description
    =====================  =============================================:
    Data                    data from ImportData()
    Freq                    frequency from dataset (must be same datatype etc)
    lin_point_shift         sometimes the linearisation is a bit weird (play around with!)
    
    return type
       Dictionary with all data you need
    
    Example:
        
        import ADS_Toolbox as ads
        
        Rename_Table= {'dBm(HB.Vout[::,1])' : "Pout_1"}
        Rename_Table['dBm(HB.Vout[::,3])'] = "Pout_3"

        # Import path
        CS_2x60 = '../Kennlinien/HB/2x60_12V_1V25.txt'

        # Import Simulation
        Data = ads.ImportData(CS_2x60, Rename=True, Rename_Table=Rename_Table, Merge=True, MergeColumn='F0') 
        
        # Calculate Harmonic balance
        HB_CS_8x60_Single = ads.Calc_HarmonicBalance_Single(D_CS_8x60, 10.0, lin_point_shift=2)      
    """   
#############################################################################   
   # Generate Frequencyset from Data
    Frequencyset = []
    for dataset in Data:
        if 'F0' in dataset:
            Frequencyset.append(dataset['F0'][0])

    # Looking for index
    Dataset_index = Frequencyset.index(Freq)
    
    # extract Dataset
    Data_Single = Data[Dataset_index]
    
    # Calculate HB
    Data_HB = []
    Data_HB.append(Data_Single)
    Data_HB = Calc_HarmonicBalance(Data_HB, debugvalue=True, lin_point_shift=lin_point_shift,
                                   lin_point_tolerance=lin_point_tolerance, forcezero=forcezero)
    
    # Extrapolate Graph
    x = np.array(Data_Single['Pin'])
    
    # calc offset
    offset_pout1 = Data_Single['Pout_1'][Data_HB['Y_Index_1'][0]] - 1*x[Data_HB['Y_Index_1'][0]]
    offset_pout2 = Data_Single['Pout_2'][Data_HB['Y_Index_2'][0]] - 2*x[Data_HB['Y_Index_2'][0]]
    offset_pout3 = Data_Single['Pout_3'][Data_HB['Y_Index_3'][0]] - 3*x[Data_HB['Y_Index_3'][0]]
        
    # Linear Functions
    Extrapolate_Pout1 = 1*x + offset_pout1
    Extrapolate_Pout2 = 2*x + offset_pout2
    Extrapolate_Pout3 = 3*x + offset_pout3

     # return type
    HarmonicBalance = {}
    HarmonicBalance['Freq'] = Freq
    HarmonicBalance['Pin'] = Data_Single['Pin']
    HarmonicBalance['Pout_1'] = Data_Single['Pout_1']
    HarmonicBalance['Extrapolation_Pout_1'] = Extrapolate_Pout1
    HarmonicBalance['Pout_2'] = Data_Single['Pout_2']
    HarmonicBalance['Extrapolation_Pout_2'] = Extrapolate_Pout2
    HarmonicBalance['Pout_3'] = Data_Single['Pout_3']
    HarmonicBalance['Extrapolation_Pout_3'] = Extrapolate_Pout3
    
    #stuff from normal HB
    HarmonicBalance['P1dB_in'] = Data_HB['P1dB_in']
    HarmonicBalance['P1dB_out'] = Data_HB['P1dB_out']
    HarmonicBalance['IP2_in'] = Data_HB['IP2_in']
    HarmonicBalance['IP2_out'] = Data_HB['IP2_out']
    HarmonicBalance['IP3_in'] = Data_HB['IP3_in']
    HarmonicBalance['IP3_out'] = Data_HB['IP3_out']
    HarmonicBalance['Psat_out'] = Data_HB['Psat_out']
     
    # return
    return HarmonicBalance

#############################################################################
###         Extract Harmonic Balance and Plot in one Diagram
#############################################################################

def Calc_HarmonicBalance(Data, debugplot=False, debugvalue=False, lin_point_shift=1, forcezero=False, lin_point_tolerance = 0.15,  **kwargs):
#############################################################################    
    """
    Calculate Harmonic Balance Values for a whole Dataset
    
    Data format should be from the function: ImportData()

    paramters              description
    =====================  =============================================:
    Data                    data from ImportData()
    debugplot               returns IP3, Psat, P1dB and Linearisation points for each dataset 
    debugvalue              returns debugvalues
    lin_point_shift         sometimes the linearisation is a bit weird (play around with!)
    
    return type
       Dictionary with Frequencies, IP3 (in + out), P1dB(in + out) and Psat_out
    
    Example:
        
        import ADS_Toolbox as ads
        
        Rename_Table= {'dBm(HB.Vout[::,1])' : "Pout_1"}
        Rename_Table['dBm(HB.Vout[::,3])'] = "Pout_3"

        # Import path
        CS_2x60 = '../Kennlinien/HB/2x60_12V_1V25.txt'

        # Import Simulation
        Data = ads.ImportData(CS_2x60, Rename=True, Rename_Table=Rename_Table, Merge=True, MergeColumn='F0') 
        
        # Calculate Harmonic balance
        HB_CS_2x60 = ads.Calc_HarmonicBalance(D_CS_2x60)      
    """    
#############################################################################    
    # Returntypes
    P1dBin = []
    P1dBout = []
    IIP2 = []
    OIP2 = []
    IIP3 = []
    OIP3 = []
    freq = []
    Psatout = []
    
    # temp: debug export
    yindex_1 = []
    yindex_2 = []
    yindex_3 = [] 
    
    # return type
    HarmonicBalance = {}

    for dataset in Data:
        
        # missing columns
        if len(dataset) < 4:
            continue
        
        # get corresponding frequency
        current_freq = dataset['F0'][0]

        # Zoom in defined frequency range
        freq_start = 0
        freq_stop = 40
        freq_Range = (freq_start < current_freq) & (current_freq < freq_stop)
 
        if not freq_Range:
            continue
        
        # get data
        if 'Pout_1' in dataset:
            y_pout1 = dataset['Pout_1']
            
        if 'Pout_2' in dataset:
            y_pout2 = dataset['Pout_2']
            
        if 'Pout_3' in dataset:
            y_pout3 = dataset['Pout_3']
            
        x = np.array(dataset['Pin'])

#############################################################################
    # try to find offset
        try:        
            if y_pout1:
                # Find Lin points for 1st Fct
                # Extract slope
                dy_pout1 = np.diff(y_pout1)/np.diff(x) 
                
                # find absolute zero crossing
                if forcezero:
                    dy_pout1 = dy_pout1-(1-lin_point_tolerance)
                
                # extract first zero
                dy_zero1 = np.where(np.diff(np.sign(dy_pout1)))[0][0]
                dy_zero1 = np.floor(dy_zero1/(lin_point_shift+1)).astype(np.int64)
                
                # x index of linearization        
                lin_xpoint1 = x[dy_zero1]
            
            if y_pout2:
                # Find Lin points for 2nd Fct
                # Extract slope
                dy_pout2 = np.diff(y_pout2)/np.diff(x) 
                
                # find absolute zero crossing
                if forcezero:
                    dy_pout2 = dy_pout2-(2-lin_point_tolerance)
                
                # extract first zero
                dy_zero2 = np.where(np.diff(np.sign(dy_pout2)))[0][0]
                dy_zero2 = np.floor(dy_zero2/(lin_point_shift+1)).astype(np.int64)
                
                # x index of linearization     
                lin_xpoint2 = x[dy_zero2] 
                
            if y_pout3:  
                # Find Lin points for 3th Fct
                # Extract slope
                dy_pout3 = np.diff(y_pout3)/np.diff(x) 
                
                # find absolute zero crossing
                if forcezero:
                    dy_pout3 = dy_pout3-(3-lin_point_tolerance)
                
                # extract first zero
                dy_zero3 = np.where(np.diff(np.sign(dy_pout3)))[0][0]
                dy_zero3 = np.floor(dy_zero3/(lin_point_shift+1)).astype(np.int64)
            
                # x index of linearization   
                lin_xpoint3 = x[dy_zero3]
            
        except:
            
            # linearize at zero
            lin_xpoint1 = 0
            lin_xpoint2 = 0
            lin_xpoint3 = 0
            
#############################################################################         
        if y_pout1:
            # y index of linearization 
            yindex1 = (np.where(x>lin_xpoint1)[0])[0] 
            # for debug export
            yindex_1.append(yindex1)
            # calc offset
            offset_pout1 = y_pout1[yindex1] - 1*x[yindex1]
            # Linear Functions
            Linear_Pout1 = 1*x + offset_pout1
                   
            # looking for Psat
            Psat_at_point = np.argmax(y_pout1)
            Psatout_local = y_pout1[Psat_at_point]
            Psatin_local = x[Psat_at_point]        
            Psatout.append(Psatout_local)

            # Calculate -1dB Compression Point nummerically
            Difference_Pout1 = np.abs(Linear_Pout1-y_pout1-1)
            Difference_Pout1_atPoint = np.argmin(Difference_Pout1)
                    
            # Search -1dB Compression Point
            P1dBin_local = x[Difference_Pout1_atPoint]
            P1dBout_local = Linear_Pout1[Difference_Pout1_atPoint]
            P1dBin.append(P1dBin_local)
            P1dBout.append(P1dBout_local)
            
            # return type
            HarmonicBalance['P1dB_in'] = P1dBin
            HarmonicBalance['P1dB_out'] = P1dBout
            HarmonicBalance['Psat_out'] = Psatout
        
 #############################################################################                  
        if y_pout2:
            # y index of linearization 
            yindex2 = (np.where(x>lin_xpoint2)[0])[0] 
            # for debug export
            yindex_2.append(yindex2)
            # calc offset
            offset_pout2 = y_pout2[yindex2] - 2*x[yindex2]
            # Linear Functions
            Linear_Pout2 = 2*x + offset_pout2
            
            # Calculate secound Order Point nummerically
            IIP2_local = (offset_pout1-offset_pout2)/1
            OIP2_local = 2*IIP2_local+offset_pout2
            IIP2.append(IIP2_local)
            OIP2.append(OIP2_local)
            
            # return type
            HarmonicBalance['IP2_in'] = IIP2
            HarmonicBalance['IP2_out'] = OIP2
        
############################################################################# 
        if y_pout3:
            # y index of linearization 
            yindex3 = (np.where(x>lin_xpoint3)[0])[0] 
            # for debug export
            yindex_3.append(yindex3)
            # calc offset
            offset_pout3 = y_pout3[yindex3] - 3*x[yindex3]
            # Linear Functions
            Linear_Pout3 = 3*x + offset_pout3

            # Calculate Thord Order Point nummerically
            IIP3_local = (offset_pout1-offset_pout3)/2
            OIP3_local = 3*IIP3_local+offset_pout3
            IIP3.append(IIP3_local)
            OIP3.append(OIP3_local)  
            
            # return type
            HarmonicBalance['IP3_in'] = IIP3
            HarmonicBalance['IP3_out'] = OIP3
             
#############################################################################             
             
        # save freq range
        freq.append(current_freq*1e9)
        HarmonicBalance['Freq'] = freq
                          
 ############################################################################# 
                
        # plot
        if debugplot:
            [fig, ax] = plt.subplots(figsize=(5, 6))
            plt.title('Frequency: %s' %current_freq)
            
            # Fundamental
            plt.plot(x, y_pout1)
            plt.plot(x, Linear_Pout1, linestyle='--', color=ax.lines[-1].get_color())
            plt.plot(x[dy_zero1+1], y_pout1[dy_zero1+1],marker='o', color=ax.lines[-1].get_color())            
            
            # 2nd Harmonic
            #plt.plot(x, y_pout2)
            #plt.plot(x, Linear_Pout2, linestyle='--', color=ax.lines[-1].get_color())
            #plt.plot(x[dy_zero2+1], y_pout2[dy_zero2+1],marker='o', color=ax.lines[-1].get_color())
            
            # 3th Harmonic
            plt.plot(x, y_pout3)
            plt.plot(x, Linear_Pout3, linestyle='--', color=ax.lines[-1].get_color())
            plt.plot(x[dy_zero3+1], y_pout3[dy_zero3+1],marker='o', color=ax.lines[-1].get_color())
            
            # Some marker
            plt.plot(P1dBin_local, P1dBout_local, marker='x', color='r')
            plt.plot(IIP2_local, OIP2_local, marker='x', color='r') 
            plt.plot(IIP3_local, OIP3_local, marker='x', color='r') 
            plt.plot(Psatin_local, Psatout_local, marker='x', color='r')  
            
            # plot properties
            #ax.set_xbound(-10,10)
            #ax.set_ybound(-40,0)
            plt.grid(True)

    if debugvalue:
        HarmonicBalance['Y_Index_1'] = yindex_1
        HarmonicBalance['Y_Index_2'] = yindex_3
        HarmonicBalance['Y_Index_3'] = yindex_3
     
    # return
    return HarmonicBalance

#############################################################################
###         Extract MSG/MAG and StabFac
#############################################################################
   
def Calc_StabGain(S_Data,  **kwargs):
#############################################################################    
    """
    Calculate Gains and Stab Fact
    
    Data format should be from converted with Conversion.py
    used also scikit-rf

    paramters              description
    =====================  =============================================:
    S_Data                  data from ImportData()
    
    return type
           K Fact, K Point (K=1), Gain, Gain dB, Mason Gain, Mason Gain in dB, Gmax, Gmax in dB and Freq
    
    Example:
        
        import Conversion as cv
        import ADS_Toolbox as ads
        
        # import s2p
        CS_2x60 = rf.Network('../Kennlinien/MAGMSG/2x60_12V_1V25.s2p')
        
        # convert
        S_CS_8x60 = cv.Network2Sparam(CS_8x60)

        # Extract Gains
        G_CS_2x60 = ads.Calc_StabGain(S_CS_2x60)      
    """   
#############################################################################   
    
    # Extract StabFact, raw data
    Det_S = S_Data['S11']*S_Data['S22'] - S_Data['S12']*S_Data['S21']
    K_Fact = 1-abs(S_Data['S11'])**2-abs(S_Data['S22'])**2+abs(Det_S)**2
    K_Fact = K_Fact / (2*(abs(S_Data['S12']*S_Data['S21'])))
    
    # Generate Stability Array (1 = MAG, -1 =MSG)
    K_Stability = np.sign(K_Fact-1)
    
    # Find index of "K_point"
    K_Point_index = np.where(np.diff(np.sign(K_Stability)))[0]

    # K_point not found
    if not K_Point_index.any():
        K_Point_index=len(K_Fact)-1
    
    # empty gain:
    Gain = np.zeros(len(Det_S))
        
    # Extract MAG/MSG
    index = 0
    for local_K in K_Stability:
        
        # MAG  (K>=1)
        if local_K == 1.0:
            MAG = abs(S_Data['S21'][index]/S_Data['S12'][index])
            MAG = MAG*(K_Fact[index]-np.sqrt(abs((K_Fact[index]**2-1))))
            Gain[index] = MAG
    
        # MSG (K<1)
        if local_K == -1.0:
            MSG = abs(S_Data['S21'][index]/S_Data['S12'][index])
            Gain[index] = MSG
            
        # iterate
        index = index + 1
    
    # Mason's Gain
    U = abs((S_Data['S21']/S_Data['S12'])-1)**2
    U = U/(2*K_Fact*abs(S_Data['S21']/S_Data['S12'])-2*np.real(S_Data['S21']/S_Data['S12']))

    # Theoretical maximum gain of the device
    Gmax = (2*U-1)+2*np.sqrt(U*(U-1))
    
    # Extract K Point
    K_Point_Freq = S_Data['Freq'][K_Point_index]
    K_Point_Gain = 10*np.log10(Gain[K_Point_index])
    
    #return types
    GainStab = {}
    GainStab['K_Fact'] = K_Fact
    GainStab['K_Stability'] = K_Stability
    GainStab['K_Point'] = [K_Point_Freq, K_Point_Gain]
    GainStab['Gain'] = Gain
    GainStab['Gain_dB'] = 10*np.log10(Gain)
    GainStab['Mason'] = U
    GainStab['Mason_dB'] = 10*np.log10(U)
    GainStab['Gmax'] = Gmax
    GainStab['Gmax_dB'] = 10*np.log10(Gmax)
    GainStab['Freq'] = S_Data['Freq']
    
    return GainStab

#############################################################################
###         Calculate Stab_Circle
#############################################################################
def Calculate_StabCircle(S_Data, freq, clip_radius = 0.0, offset_circle=0,  **kwargs):
#############################################################################    
    """
    Calculate Stability Circle at frequency
    
    Data format should be from converted with Conversion.py
    used also scikit-rf

    paramters              description
    =====================  =============================================:
    S_Data                  data from ImportData()
    freq                    selected frequency in Hz
    clip_radius             (option) cut circle outside smith chart
    offset_circle           (option) shift circle start point in degree
    
    return type
           Freq, Output Circle, Input Circle
    
    Example:
        
        import Conversion as cv
        import ADS_Toolbox as ads
        
        # import s2p
        CS_2x60 = rf.Network('../Kennlinien/MAGMSG/2x60_12V_1V25.s2p')
        
        # convert
        S_CS_8x60 = cv.Network2Sparam(CS_8x60)

        # Extract Gains
        freq =15e9 #in Hz
       
        # find circles
        Circles = ads.Calculate_StabCircle(S_CS_8x60, freq)      
    """   
#############################################################################       
    # Frequency index
    freq_index = np.argmin(np.abs(S_Data['Freq']-freq))
    
    # Extract Paramters
    S11 = S_Data['S11'][freq_index]
    S12 = S_Data['S12'][freq_index]
    S21 = S_Data['S21'][freq_index]
    S22 = S_Data['S22'][freq_index]
    
    # Determinante of the Twoport
    det_S = S11*S22 - S12*S21
    
    # full circle
    offset_circle = offset_circle*(2*np.pi)/360
    start_circle = np.pi+offset_circle
    stop_circle = -np.pi+offset_circle
    t_full = np.linspace(start_circle,stop_circle,360)
    
############################################################################# 
    # Calculate Output Circle
    output_circle_radius = np.abs((S12*S21)/(np.abs(S22)**2 - np.abs(det_S)**2))
    output_circle_center = np.conj(S22-(det_S*np.conj(S11)))/(np.abs(S22)**2- np.abs(det_S)**2)
    
    output_circle = []
    
    if clip_radius:
        # calculate clipping box
        for t in t_full:
            mag_circle = np.abs(output_circle_radius*np.exp(1j*t) + output_circle_center)
        
            if  mag_circle < clip_radius:
                output_circle.append(output_circle_radius*np.exp(1j*t) + output_circle_center)
    else:
        
        output_circle = output_circle_radius*np.exp(1j*t_full) + output_circle_center
        
############################################################################# 
    # Calculate Input Circle
    input_circle_radius = np.abs((S12*S21)/(np.abs(S11)**2 - np.abs(det_S)**2))
    input_circle_center = np.conj(S11-(det_S*np.conj(S22)))/(np.abs(S11)**2- np.abs(det_S)**2)
    
    input_circle = []
    
    if clip_radius:
        # calculate clipping box
        for t in t_full:
            mag_circle = np.abs(input_circle_radius*np.exp(1j*t) + input_circle_center)
        
            if  mag_circle < clip_radius:
                input_circle.append(input_circle_radius*np.exp(1j*t) + input_circle_center)
    else:
        
        input_circle = input_circle_radius*np.exp(1j*t_full) + input_circle_center
        
############################################################################# 
    # add Dataset
    Data = {}
    
    # Add Frequency, Standart and MixedMode row
    Data['freq'] = S_Data['Freq'][freq_index]  
    Data['Output_Circle'] = output_circle  
    Data['Input_Circle'] = input_circle 
         
    return Data         
        
#############################################################################
###         Calculate Mixed Mode Sparam
#############################################################################
def MixedModeSparam(Sim_file, elementvar='var("S")', complexdata=False,  ZParam=False, Z0=50, **kwargs): 
#############################################################################    
    """
    Calculate MixedMode Parameters
    
    Data format should be from converted with Conversion.py
    used also scikit-rf

    paramters              description
    =====================  =============================================:
    Sim_file                Import Data via ImportData-Fct
    elementvar             (option) name of Sparam-Elements
    complexdata            (option) import complex data
    ZParam                 (option) calculate Z Parameters
    Z0                     (option) charateristic impedance   
    
    return type
           freq, Standart, MixedMode, MixedMode_Phase and MixedMode_dB Matrix 
    
    Example:
        
        import ADS_Toolbox as ads
        
        Sim_file ='../../Schaltungen/ADF102GH020_CHO/ADF102GH020_CHO_differential_with_L.txt'

        Simulation = ads.MixedModeSparam(Sim_file, complexdata=True)

        # .....
        # Figure etc for plot
        
        ax1.semilogx(Simulation['freq'], Simulation['MixedMode_dB']['Sdd21'])
    
    """   
#############################################################################    
    Matrix_4Port = ImportData(Sim_file, Merge=True, MergeColumn='freq', complexdata=complexdata)
    Matrix_4Port = Matrix_4Port[0]
        
    # Mixed Mode Matrix
    Mixed_Matrix = [[1, -1, 0, 0], [0, 0, 1, -1], [1, 1, 0, 0], [0, 0, 1, 1]]
    Mixed_Matrix = (1/np.sqrt(2)) * np.array(Mixed_Matrix)
    Mixed_Matrix_invers = np.linalg.inv(np.matrix(Mixed_Matrix))
    
        
    # frequency array
    index_freq = np.linspace(0,len(Matrix_4Port['freq'])-1, num=len(Matrix_4Port['freq']))
    index_freq = index_freq.astype(int)
        
    # Generate Standart Matrix
    index_column = [1,2,3,4]
    index_row =  [1,2,3,4]
    Sparam_Standart = np.zeros((len(index_row),(len(index_column))), dtype=np.complex)
    
    # generate empty Standart and Mixed Mode
    Standart = {}
    MixedMode = {}
    MixedMode_dB = {}
    MixedMode_Phase = {}
    
    # Z param calc
    if ZParam:
        ZMatrix = {}
        ZMatrix_mag = {}
    
    for row in index_row:
        for col in index_column:
            Standart['S' + str(col) + str(row)] = []
            MixedMode['S' + str(col) + str(row)] = [] 
            MixedMode_dB['S' + str(col) + str(row)] = []   
            MixedMode_Phase['S' + str(col) + str(row)] = []   
            
            # Z param calc
            if ZParam:
                ZMatrix['Z' + str(col) + str(row)] = []
                ZMatrix_mag['Z' + str(col) + str(row)] = []
                    
    # iterate Frequency
    for freq in index_freq:
                
        # Iterate Matrix
        for row in index_row:
            for col in index_column:
                # Import to Matrix
                Element = Matrix_4Port[elementvar + '('+ str(col) +','+ str(row) +')('+ str(col) +','+ str(row) +')']
                Sparam_Standart[col-1][row-1] = Element[freq]
         
        # Calculate Mixed Mode Matrix
        Sparam_MixedMode = np.dot(Sparam_Standart, Mixed_Matrix_invers)
        Sparam_MixedMode =  np.dot(Mixed_Matrix, Sparam_MixedMode)
        Sparam_MixedMode = np.asarray(Sparam_MixedMode)
        
        # Calculate Z Matrix
        if ZParam:
            
            # generate some useful arrays
            M_ident = np.identity(4)
            M_sqrtZ0 = np.identity(4)*np.sqrt(Z0)
            
            # Generate ZMatrix
            Z_Invers = np.linalg.inv(np.matrix(M_ident-Sparam_MixedMode))
            Z_Normal = np.matrix((M_ident+Sparam_MixedMode))
            
            Zparam_Matrix = np.dot(Z_Invers,M_sqrtZ0)
            Zparam_Matrix = np.dot(Z_Normal,Zparam_Matrix)
            Zparam_Matrix = np.dot(M_sqrtZ0,Zparam_Matrix)
            Zparam_Matrix = np.asarray(Zparam_Matrix)
            
            
        # fill Std and MixedMode Matrix   
        for row in index_row:
            for col in index_column:
                Standart['S' + str(col) + str(row)].append(Sparam_Standart[col-1][row-1])
                MixedMode['S' + str(col) + str(row)].append(Sparam_MixedMode[col-1][row-1])
                MixedMode_dB['S' + str(col) + str(row)].append(20*np.log10(np.absolute(Sparam_MixedMode[col-1][row-1])))
                MixedMode_Phase['S' + str(col) + str(row)].append(np.angle(Sparam_MixedMode[col-1][row-1], deg=True))
                
                # Z matrix
                if ZParam:
                    ZMatrix['Z' + str(col) + str(row)].append(Zparam_Matrix[col-1][row-1])
                    ZMatrix_mag['Z' + str(col) + str(row)].append(np.abs(Zparam_Matrix[col-1][row-1]) )                  
      
            
    # add Dataset
    Data = {}
    
    # Add Frequency, Standart and MixedMode row
    Data['freq'] = Matrix_4Port['freq']   
    Data['Standart'] = Standart  
    Data['MixedMode'] = MixedMode 
    Data['MixedMode_dB'] = MixedMode_dB
    Data['MixedMode_Phase'] = MixedMode_Phase
    
    if ZParam:
        Data['ZMatrix'] = ZMatrix 
        Data['ZMatrix_mag'] = ZMatrix_mag
         
    return Data 

#############################################################################
###         Import S2P-File
#############################################################################
def ImportS2P(File, ZParam = False, HParam=False, PhaseOffs=True,  **kwargs):
#############################################################################    
    """
    Import S2P File with different parameters
    
    Data format should be from converted with Conversion.py
    used also scikit-rf

    paramters              description
    =====================  =============================================:
    S_Data                  data from ImportData()
    freq                    selected frequency in Hz
    clip_radius             (option) cut circle outside smith chart
    
    return type
           Freq, Output Circle, Input Circle
    
    Example:
        
        import Conversion as cv
        import ADS_Toolbox as ads
        
        # import s2p
        CS_2x60 = rf.Network('../Kennlinien/MAGMSG/2x60_12V_1V25.s2p')
        
        # convert
        S_CS_8x60 = cv.Network2Sparam(CS_8x60)

        # Extract Gains
        freq =15e9 #in Hz
       
        # find circles
        Circles = ads.Calculate_StabCircle(S_CS_8x60, freq)      
    """   
#############################################################################      
    # import Network
    Network = rf.Network(File)
    
    # Extract complex S Param
    Complex_Network = {}
    Complex_Network['S11'] = Network.s11.s.reshape(-2)
    Complex_Network['S12'] = Network.s12.s.reshape(-2)
    Complex_Network['S21'] = Network.s21.s.reshape(-2)
    Complex_Network['S22'] = Network.s22.s.reshape(-2)
    
    # Extract Z-Param:
    if ZParam:
        Complex_Network['Z11'] = Network.z[:,0,0]
        Complex_Network['Z12'] = Network.z[:,0,1]
        Complex_Network['Z21'] = Network.z[:,1,0]
        Complex_Network['Z22'] = Network.z[:,1,1]
    
    # Extract H-Param:
    if HParam:
        
        # Equations from HHHS Building Blocks - Active Devices Page 29
        det_S = Complex_Network['S11']*Complex_Network['S22']
        det_S = det_S-Complex_Network['S12']*Complex_Network['S21']
        
        denom_H = -det_S - Complex_Network['S11'] + Complex_Network['S22'] + 1
    
        # H Param
        Complex_Network['H11'] = (det_S + Complex_Network['S11'] + Complex_Network['S22'] + 1)
        Complex_Network['H11'] = Complex_Network['H11']/denom_H
        
        Complex_Network['H12'] = (2*Complex_Network['S12'])/denom_H
                       
        Complex_Network['H21'] = (-2*Complex_Network['S21'])/denom_H
                       
        Complex_Network['H22'] = (det_S - Complex_Network['S11'] + Complex_Network['S22']  + 1)
        Complex_Network['H22'] = Complex_Network['H22']/denom_H        
    
    
    # Generate Return Type
    Data = {}
    
    for Cleantag in Complex_Network:
        Data[Cleantag] = Complex_Network[Cleantag]
        Data[Cleantag + "_dB"] = 20*np.log10(np.absolute(Complex_Network[Cleantag]))
        
        # extract phase in degree and unwrapped
        phase = np.angle(Complex_Network[Cleantag])
        phase = np.unwrap(phase)
        
        # determine phase offset
        phase_offset = 0
        angle_offset = 0
        
        if PhaseOffs:
            phase_offset = phase[0]
            angle_offset = np.angle(Complex_Network[Cleantag][0])
            
        Data[Cleantag + "_Phase"] = (phase-phase_offset)*360/(2*np.pi)
        Data[Cleantag + "_Angle"] = np.angle(Complex_Network[Cleantag])-angle_offset  

    # Environment
    Data['Z0'] = Network.z0  
    Data['freq'] = Network.frequency.f


    return Data