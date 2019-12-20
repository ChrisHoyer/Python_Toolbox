# Measure Sampling Scope
# C.Hoyer (christian.hoyer1@tu-dresden.de)
# Equipment: RS_FSU67, RS_FSW67

###############################################################
import numpy as np
import visa, time, datetime, csv, logging, os

#############################################################################
#           Import used Instruments
#############################################################################
import _VISAInstruments.Instrument
from  _VISAInstruments.Tektronix_DSA8300 import Tektronix_DSA8300

#############################################################################
#           Generate Output Format
#############################################################################
def Generate2CSV (data, setting, filename):
    
    # Generate File
    file = open(filename,'w')
    
    # Write Settings
    file. write(setting + '\n')
    
    # Write Data
    file. write(data)
    
    # close file
    file.close()

#############################################################################
#          Init
#############################################################################

# Instrument definition
dsa = Tektronix_DSA8300(addr= "GPIB1::1::INSTR") #DSO

# Set Trigger Settings
dsa.set_trigger_source("EXTdirect")
dsa.set_trigger_slope("FALL")
dsa.set_trigger_level(200e-3)

# Horizontal Axis, Main Timebase
dsa.set_MainTB_scale(20e-12)

# Show channels 1 and 2
dsa.set_Channel_select(1, "ON")
dsa.set_Channel_select(2, "ON")

#############################################################################
#          Export Current 
#############################################################################

# Export Data
dsa.set_Waveform_DataFormat("ASCIi")

# stop 
dsa.set_acq_state("STOP")

dsa.set_Waveform_Source_TB("CH1", "MAIn")
data_CH1 = dsa.get_Waveform_Curve()
setting_CH1 = dsa.get_Waveform_Preamble()

dsa.set_Waveform_Source_TB("CH2", "MAIn")
data_CH2 = dsa.get_Waveform_Curve()
setting_CH2 = dsa.get_Waveform_Preamble()

# start 
dsa.set_acq_state("RUN")

# Save Data
Generate2CSV(data_CH1, setting_CH1, "2019_12_13_SamplingScope_Coupled_CH1.csv")
Generate2CSV(data_CH2, setting_CH2, "2019_12_13_SamplingScope_Coupled_CH2.csv")

print("Exported: Single Run")

#############################################################################
#          Export Database with long-term
#############################################################################

# Set Waveform Database
dsa.set_WfmDB(1, "CH1")
dsa.set_WfmDB(2, "CH2")

# Show Waveform Database
dsa.show_WfmDB(1, "ON")
dsa.show_WfmDB(2, "ON")

# Set Length
dsa.set_WfmDB_persis(1, 500, "VARPersist")
dsa.set_WfmDB_persis(2, 500, "VARPersist")

