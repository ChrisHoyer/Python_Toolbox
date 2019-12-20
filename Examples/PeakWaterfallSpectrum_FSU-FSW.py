# Measure Spectrum over Time
# C.Hoyer (christian.hoyer1@tu-dresden.de)
# Equipment: RS_FSU67, RS_FSW67

###############################################################
import numpy as np
import visa, time, datetime, csv, logging, os

#############################################################################
#           Import used Instruments
#############################################################################
import _VISAInstruments.Instrument
from  _VISAInstruments.RS_FSW_FSU import RS_FSW_FSU


#############################################################################
#           Measurement Functions
#############################################################################
def init_SpectrumAnalyser(sa):
    #sa.reset()
    
    # sweep mode
    sa.set_mode_CW()
    
    # Frequency Spectrum
    sa.set_freq_center(24.6165e9)
    sa.set_freq_span(10e6)
    sa.set_resolution_bw(100e3)     
    sa.set_sweep_points(1001)
    sa.set_unit("DBM")
    
    # Sweep
    sa.set_sweep_time(50e-3)
    
    # Marker Functions
    sa.set_marker_count()
    
#############################################################################
def measure_SinglePeak(sa, starttime):
          
    # find maximum
    sa.move_marker_peak(1)
     
    ydata = sa.get_marker_y(1)
    xdata = sa.get_marker_x(1)     
     
    # generate current time stamp
    meas_time = time.time() - starttime
     
    # return data
    return [meas_time, xdata, ydata]

#############################################################################
#          	Main Loop
#############################################################################    
   
# Instrument definition
fsu = RS_FSW_FSU(addr= "TCPIP0::172.31.228.61::inst0::INSTR") # FSU
fsw = RS_FSW_FSU(addr= "TCPIP0::172.31.228.85::inst0::INSTR") # FSW

# Instrument Init
init_SpectrumAnalyser(fsu)
init_SpectrumAnalyser(fsw)   

# Some Settling Time
time.sleep(5)

# Perform N Samples
N = 10

# Write FSU Data
file_fsu = open("F:\2019_12_12_Coupled_LongTime_FSU.csv",'w')
file_fsu.write("Time[s]" + ", Freq[Hz]" + ", Mag[dBm]" + '\n')

# Write FSW Data
file_fsw = open("F:\2019_12_12_Coupled_LongTime_FSW.csv",'w')
file_fsw.write("Time[s]" + ", Freq[Hz]" + ", Mag[dBm]" + '\n')

#print(xData)
print(yData)

fsw.set_mode_CW()

# Start of Measurement
print("Starting Measurement")

# Start Time
start_time = time.time()

# Measurement Loop
for i in range(N):
    
    # Measurement
    data_fsu = measure_SinglePeak(fsu, start_time)
    data_fsw = measure_SinglePeak(fsw, start_time)
    
    
    # Write Data in File
    file_fsu.write(str(data_fsu).lstrip('[').rstrip(']') + '\n')
    file_fsw.write(str(data_fsw).lstrip('[').rstrip(']') + '\n')
     
    # Print Data
    #print('FSU: ' + str(data_fsu[0]) + ',' + str(data_fsu[1]/1e9) + 'GHz ,' + str(data_fsu[2]) + 'dBm' )
    #print('FSW: ' + str(data_fsw[0]) + ',' + str(data_fsw[1]/1e9) + 'GHz ,' + str(data_fsw[2]) + 'dBm' )
    
    # Wait inbetween
    #time.sleep(0.5)

file_fsu.close() 
file_fsw.close() 

# End of Measurement
print("Finished Measurement")



    


    
    
    
    


