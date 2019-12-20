# Measure Spectrum over Time
# C.Hoyer (christian.hoyer1@tu-dresden.de)
# Equipment: RS_FSU67, RS_FSW67, DSA8300

###############################################################
import numpy as np
import visa, time, datetime, csv, logging, os

#############################################################################
#           Import used Instruments
#############################################################################
import sys
sys.path.insert(1, '//iee.local/CCN/Homes/christian.hoyer/Documents/Projekte/_Templates/Python/_Labinstruments')

import Instrument
from  Instruments.RS_FSW_FSU import RS_FSW_FSU
from  Instruments.Tektronix_DSA8300 import Tektronix_DSA8300

#############################################################################
#           Measurement Functions
#############################################################################
def init_SpectrumAnalyser(sa):
    #sa.reset()
    
    # sweep mode
    sa.set_mode_CW()
    
    # Frequency Spectrum
    sa.set_freq_center(24.47e9)
    sa.set_freq_span(100e6)
    sa.set_resolution_bw(100e3)     
    sa.set_sweep_points(1001)
    sa.set_unit("DBM")
	
	# reference level
    sa.set_RefLevel(0)
    
    # Sweep
    sa.set_sweep_time(50e-3)
    
    # Marker Functions
    sa.set_marker_count()
    
#############################################################################
def init_SamplingScope(sa):

	# Set Trigger Settings
	sa.set_trigger_source("EXTdirect")
	sa.set_trigger_slope("FALL")
	sa.set_trigger_level(200e-3)

	# Horizontal Axis, Main Timebase
	sa.set_MainTB_scale(20e-12)

	# Show channels 1 and 2
	sa.set_Channel_select(1, "ON")
	sa.set_Channel_select(2, "ON")
	
	# Export Data
	dsa.set_Waveform_DataFormat("ASCIi")
	
	# Start Device
	dsa.set_acq_state("RUN")
    
#############################################################################
def measure_SinglePeak(sa, starttime):
    
    # Generate Output File
    Return_Dict = {}
          
    # find maximum
    sa.move_marker_peak(1)
     
    Return_Dict["yData"] = sa.get_marker_y(1)
    Return_Dict["xData"] = sa.get_marker_x(1)     
     
    # generate current time stamp
    Return_Dict["time"] = time.time() - starttime
     
    # return data
    return Return_Dict
	
#############################################################################
def measure_SamplingScope(sa, ch, starttime):
    
    # Generate Output File
    Return_Dict = {}
   
    # get CH Data
    dsa.set_Waveform_Source_TB(ch, "MAIn")
    
    # sub measurement start-time
    Return_Dict["Measuring_Time"] = time.time() - starttime
    
    # get data
    data_raw = dsa.get_Waveform_Curve()
	
	# generate data list
    data = []
    data.append([float(s) for s in data_raw.split(',')])
    
    Return_Dict["Data"] = data
	
	# Channel Settings
    Return_Dict["Setting"] = dsa.get_Waveform_Preamble()
    
    # return data
    return Return_Dict
	
#############################################################################
#           Measurement Loop
############################################################################
def measure(sa1, sa2, dsa, start_time):

    # Generate Output File
    Return_Dict = {}
    
    # Measure Frequency Peaks
    Return_Dict["Frequency_1"] = measure_SinglePeak(sa1, start_time)
    Return_Dict["Frequency_2"] = measure_SinglePeak(sa2, start_time)

    # Measure Sampling Scope
    dsa.set_acq_state("STOP")
    Return_Dict["Sampling_CH1"] = measure_SamplingScope(dsa, "CH1", start_time)
    Return_Dict["Sampling_CH2"] = measure_SamplingScope(dsa, "CH2", start_time)
    dsa.set_acq_state("RUN")	
    
    return Return_Dict

#############################################################################
#           Measurement Main
#############################################################################

# Instrument definition
fsu = RS_FSW_FSU(addr= "TCPIP0::172.31.228.61::inst0::INSTR") # FSU
fsw = RS_FSW_FSU(addr= "TCPIP0::172.31.228.85::inst0::INSTR") # FSW
dsa = Tektronix_DSA8300(addr= "GPIB1::1::INSTR") #DSO

# Instrument Init
init_SpectrumAnalyser(fsu)
init_SpectrumAnalyser(fsw) 
init_SamplingScope(dsa)

# Settling Time
time.sleep(1)

# Begin Measurement
print("First Sampling Points:")
data_fsu = measure_SinglePeak(fsu, time.time())
data_fsw = measure_SinglePeak(fsw, time.time())
# Print Data
print('FSU: ' + str(data_fsu["xData"]/1e9) + 'GHz   ,' + str(data_fsu["yData"]) + 'dBm' )
print('FSW: ' + str(data_fsw["xData"]/1e9) + 'GHz   ,' + str(data_fsw["yData"]) + 'dBm' )

#############################################################################

# cycles
N = 100
filename = "2019_12_13_Sampling_Frequency_ShortTime_Length" + str(N)

# Generate Output Data
ExportFile = open(filename + ".csv",'w')
Numpy_Export = {}

# Start Time
start_time = time.time()

#############################################################################

# Measurement Loop
for i in range(N):
    
    # Do Measurement
    data = measure(fsu, fsw, dsa, start_time)
    
    # write to file
    ExportFile.write(str(data) + '\n')
    Numpy_Export[i] = data


# Generate File
np.save(filename, Numpy_Export)
ExportFile.close()