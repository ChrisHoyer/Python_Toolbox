# packages
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import csv

# file for plotting
csv_file = "//iee.local/CCN/Homes/christian.hoyer/Downloads/Freq_Volt.csv"
delimiter = ',' 


# Import SWP Matrix
csv_reader = csv.reader(open(csv_file),  delimiter=delimiter)
    
# generate new Block
line_count = 0
Block = {}
             
# Import Block
for row in csv_reader:
            
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