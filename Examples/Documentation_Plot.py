
import sys
sys.path.insert(1, '//iee.local/CCN/Homes/christian.hoyer/Documents/Projekte/_Templates/Python')
       
import PCB_Toolbox as pcb 
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.patches as patches
import numpy as np

from scipy.optimize import curve_fit

        
# import csv
#file = 'Outline Vertices.csv'     
# calculate
#Area = pcb.CSV2Area(file, draw_polygon=True)  
#print('This Polygon has an Area of ' + str(Area) + ' mm2 or ' + str(Area/100) + ' cm2')    

#############################################################################
PreBias_Frequency = 24

# Plot VCO Transfer Function
FVCO = pcb.CSV2Dict("Freq_Volt.csv", delimiter=",")
FVCO["TuningLaw"] = np.asarray(FVCO["TuningLaw"])/1e9
FVCO["TuningLaw"] = FVCO["TuningLaw"].tolist()

# Plot VCO Transfer Function
KVCO = pcb.CSV2Dict("Kv_Volt.csv", delimiter=",")
KVCO["TuningSensitivity"] = np.asarray(KVCO["TuningSensitivity"])/1e6
KVCO["TuningSensitivity"] = KVCO["TuningSensitivity"].tolist()

# Find Point
FVCO_PreV = pcb.FindPoint_FitFct(np.asarray(FVCO["TuningLaw"]), 
                                       np.asarray(FVCO["Voltage"]), 
                                       PreBias_Frequency, 5)


KVCO_Sens = pcb.FindPoint_FitFct(np.asarray(KVCO["Voltage"]), 
                                 np.asarray(KVCO["TuningSensitivity"]), 
                                 FVCO_PreV, 5)

# Generate Marker
FVCO_PreBias = [[FVCO["Voltage"][0],FVCO_PreV, FVCO_PreV],
                [PreBias_Frequency, PreBias_Frequency, np.min(FVCO["TuningLaw"])]]


KVCO_PreBias = [[KVCO["Voltage"][0],FVCO_PreV, FVCO_PreV],
                [KVCO_Sens, KVCO_Sens, np.min(KVCO["TuningSensitivity"])]]


#############################################################################
# Calculate Linear Response
omega_0 = PreBias_Frequency*1e9 - KVCO_Sens*1e6*1*FVCO_PreV - KVCO_Sens*0.85*1*1*(1.6/2)
omega_0 = omega_0

omega_k = omega_0 + KVCO_Sens*1e6*1*np.asarray(FVCO["Voltage"]) + KVCO_Sens*0.85*1*1*(1.6/2)
omega_k = omega_k/1e9
#############################################################################

# Plot FVCO
f,ax = plt.subplots(figsize=(5, 5))

#plt.plot(KVCO["Voltage"], KVCO["TuningSensitivity"])
plt.plot(FVCO["Voltage"], FVCO["TuningLaw"], 'o-', linewidth=2, markersize=5, label="specified dataset")
plt.plot(FVCO_PreBias[0], FVCO_PreBias[1], 'k--')

# plot inst freq
plt.plot(FVCO["Voltage"], omega_k, ':', linewidth=2, label=r'linearized $\omega_k$')


# limits
ax.set_xlim(np.min(FVCO["Voltage"]),np.max(FVCO["Voltage"]))
ax.set_ylim(np.min(FVCO["TuningLaw"]),np.max(FVCO["TuningLaw"]))

# label
ax.set_ylabel('$f_{VCO}$')
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g GHz'))
ax.set_xlabel('$V_{tune}$')
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g V'))

# text
ax.plot(FVCO_PreV,PreBias_Frequency,'Hr', markersize=7)
ax.text(FVCO_PreV+0.5, PreBias_Frequency-0.5, 'pre-bias for %.2f GHz at %.2f V' 
        %(PreBias_Frequency, FVCO_PreV),
        bbox=dict(color='white', alpha=0.75))

ax.legend()

# box 
rect = patches.Rectangle((0,20),4,10,facecolor='#e6e6e6')
ax.add_patch(rect)
    
ax.annotate("", xy=(0, 27), xytext=(4, 27), arrowprops=dict(arrowstyle="<-",linewidth=2))
ax.text(0.5, 27.2, '$\mathrm{V}_\mathrm{tune,max}$', fontsize=15)
    
# grid
plt.minorticks_on()
plt.grid(True, which='major')
ax.grid(which='minor', alpha=1, linestyle=':', linewidth=1)
ax.grid(which='major', alpha=1, linewidth=1.2)
plt.show

#plt.savefig("VCO_Freq_Voltage.pdf", dpi=120, bbox_inches='tight')

#############################################################################

# Plot FVCO
f,ax = plt.subplots(figsize=(5, 5))

#plt.plot(KVCO["Voltage"], KVCO["TuningSensitivity"])
plt.plot(KVCO["Voltage"], KVCO["TuningSensitivity"], 'o-', linewidth=2, markersize=5)
plt.plot(KVCO_PreBias[0], KVCO_PreBias[1], 'k--')

# limits
ax.set_xlim(np.min(KVCO["Voltage"]),np.max(KVCO["Voltage"]))
ax.set_ylim(np.min(KVCO["TuningSensitivity"]),np.max(KVCO["TuningSensitivity"]))

# label
ax.set_ylabel('$K_{VCO}$')
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g MHz/V'))
ax.set_xlabel('$V_{tune}$')
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g V'))

# text
ax.plot(FVCO_PreV,KVCO_Sens,'Hr', markersize=7)
ax.text(FVCO_PreV+0.5, KVCO_Sens+100, 'sensitivity of %s MHz/V at pre-bias'
        %(np.around(KVCO_Sens, decimals=2)),
         bbox=dict(color='white', alpha=0.75))

# box 
rect = patches.Rectangle((0,200),4,2000,facecolor='#e6e6e6')
ax.add_patch(rect)
    
ax.annotate("", xy=(0, 1400), xytext=(4, 1400), arrowprops=dict(arrowstyle="<-",linewidth=2))
ax.text(0.75, 1450, '$\mathrm{V}_\mathrm{tune,max}$', fontsize=15)

# grid
plt.minorticks_on()
plt.grid(True, which='major')
ax.grid(which='minor', alpha=1, linestyle=':', linewidth=1)
ax.grid(which='major', alpha=1, linewidth=1.2)
plt.show

#plt.savefig("VCO_Sensitivity_Voltage.pdf", dpi=120, bbox_inches='tight')
