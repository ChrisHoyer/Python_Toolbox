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
    sa.set_mode_Sweep()
    
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
#          	Main Loop
#############################################################################    
   
# Instrument definition
fsu = RS_FSW_FSU(addr= "TCPIP0::172.31.228.61::inst0::INSTR") # FSU

# Instrument Init
init_SpectrumAnalyser(fsu)  




# End of Measurement
print("Finished Measurement")



    


    
    
    
    


