import sys
sys.path.insert(1, '//iee.local/CCN/Homes/christian.hoyer/Documents/Projekte/_Templates/Python')
import PCB_Toolbox as pcb 
import numpy as np

# import csv-file
file = 'SPI Bus 2019_10_30/scope_1.csv'
file = 'SPI Bus 2019_10_30/WA000001.CSV'
file = 'Startup_Procedure.CSV'
Data = pcb.CSV2Dict(file, delimiter=',',headerline_ends=10)

# Generate Stream
stream = pcb.Digitalize_Data(Data['CH3'],Data['CH1'], plot=True, edge_trigger='rising')

# neglect first item
stream = np.delete(stream,0)

# generate 8-bit length values
n = 16
Byte_Stream = [stream[i * n:(i + 1) * n] for i in range((len(stream) + n - 1) // n )]  

for Transmission_Block in Byte_Stream:
    n = 4
    Short = [Transmission_Block[i * n:(i + 1) * n] for i in range((len(Transmission_Block) + n - 1) // n )]  
    print(str(Short) + "\n")