
import sys
sys.path.insert(1, '//iee.local/CCN/Homes/christian.hoyer/Documents/Projekte/_Templates/Python')
       
import PCB_Toolbox as pcb 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np



if False:
    # Plot VCO Transfer Function
    XOR_SE = pcb.CSV2Dict("XOR_to_SingleEnded_1GHz.txt", delimiter="\t")
    
    # limit time
    XOR_SE["time"] = XOR_SE["time"][0:int(np.size(XOR_SE["time"])/3)]
    Label = XOR_SE["time"][int(np.size(XOR_SE["time"])/16)]
    
#############################################################################
    # Plot FVCO
    
    X_label = ['', 's']
    
    Ylabel_1 = [r'$\mathrm{V}_\mathrm{PD,in}$', 'V']
    Plot_1 = [[XOR_SE["time"], XOR_SE["V(tau_fb)"], "Feedback"],
               [XOR_SE["time"], XOR_SE["V(tau_1)"], "Cross Coupling"]]
    
    Ylabel_2 = [r'$\mathrm{V}_\mathrm{PD,out}$', 'V']
    Plot_2 = [[XOR_SE["time"], XOR_SE["V(a)"], r'$\mathrm{A}$'],
               [XOR_SE["time"], XOR_SE["V(a_inv)"], r'$\overline{\mathrm{A}}$']]
    
    Ylabel_3 = [r'$\overline{\mathrm{V}_\mathrm{\mathrm{A}-\overline{\mathrm{A}}}}$', 'V']
    Plot_3 = [[XOR_SE["time"], XOR_SE["V(a_se)"], "Single Ended"]]
    
    # Generate Plot
    plt.figure(figsize=(7.5,12))
    
    ax1 = plt.subplot(311)
    pcb.Linear_Plot(ax1, Plot_1, X_label, Ylabel_1)
    
    ax2 = plt.subplot(312)
    pcb.Linear_Plot(ax2, Plot_2, X_label, Ylabel_2)
    
    ax3 = plt.subplot(313)
    pcb.Linear_Plot(ax3, Plot_3, X_label, Ylabel_3)
    ax3.text(Label, 1, 'average voltage %.2f V' 
             %np.mean(XOR_SE["V(a_se)"]),
               bbox=dict(color='white', alpha=0.85))
    
    plt.show
    
    #plt.savefig("XOR_to_SingleEnded_1GHz.pdf", dpi=120, bbox_inches='tight')
    
#############################################################################
    
if False:
    # Plot VCO Transfer Function
    XOR_AC = pcb.CSV2Dict("XOR_to_SingleEnded_AC.txt", delimiter="\t", complexdelimiter=",")
    XOR_Magnitude = np.abs(XOR_AC["V(a_se_att)"])
    # minus in magnitue!
    Phase = np.unwrap(np.angle(XOR_AC["V(a_se_att)"]))-180
    
    XOR_AC_470 = pcb.CSV2Dict("XOR_to_SingleEnded_AC_432R.txt", delimiter="\t", complexdelimiter=",")
    XOR_Magnitude_470 = np.abs(XOR_AC_470["V(a_se)"])
    # minus in magnitue!
    Phase_470 = np.unwrap(np.angle(XOR_AC_470["V(a_se)"]))-180   
    
    # 3dB BW
    #3dB = XOR_SE["time"][int(np.size(XOR_SE["time"])/16)]
    
#############################################################################
    # Plot FVCO
    
    X_label = ['', 'Hz']
    
    Ylabel_1 = [r'$|~\overline{\mathrm{V}_\mathrm{\mathrm{A}-\overline{\mathrm{A}}}}~|$', 'V']
    Plot_1 = [[XOR_AC["Freq."], XOR_Magnitude, r'Magnitue $50~\Omega$']]
    
    Plot_11 = [[XOR_AC_470["Freq."], XOR_Magnitude_470, r'Magnitue $423~\Omega$']]
    
    Ylabel_2 = [r'$\angle~\overline{\mathrm{V}_\mathrm{\mathrm{A}-\overline{\mathrm{A}}}}~$', '$^\circ$']
    
    Plot_2 = [[XOR_AC["Freq."], Phase, r'Phase $50~\Omega$']]    
    Plot_22 = [[XOR_AC_470["Freq."], Phase_470, r'Phase $423~\Omega$']]   
    
    # Generate Plot
    plt.figure(figsize=(12,5))
    
    ax1 = plt.subplot(111)
    pcb.SemiLogX_Plot(ax1, Plot_1, X_label, Ylabel_1, Ylim=[0,2.5])
    pcb.SemiLogX_Plot(ax1, Plot_11, X_label, Ylabel_1, Ylim=[0,2.5])
    
    ax2 = ax1.twinx() 
    pcb.SemiLogX_Plot(ax2, Plot_2, X_label, Ylabel_2, TwinX=ax1, linestyle='--', color='C0')    
    pcb.SemiLogX_Plot(ax2, Plot_22, X_label, Ylabel_2, TwinX=ax1, linestyle='--', color='C1')    
    
    # align axis
    plt.show
    
    plt.savefig("XOR_to_SingleEnded_AC.pdf", dpi=120, bbox_inches='tight')
  

#############################################################################
    
if False:
    # Plot VCO Transfer Function
    JFET = pcb.CSV2Dict("JFET_Characteristic.txt", delimiter="\t", blockheader=False, dictkey=True)
    JFET_lin = pcb.CSV2Dict("JFET_Characteristic_lin.txt", delimiter="\t", blockheader=False, dictkey=True)    
    
    # pop empty set
    JFET.pop('')
    JFET_lin.pop('')
    
    # Plot list
    Plot = []
    Plot_lin = []
    Plot1 = []
    Plot1_lin = []
    
    # rename dict keys
    for dict_key in JFET.keys():
        
        # generate list
        if(len(dict_key.split(' ')) > 3):

            # generate list
            Plot.append([JFET[dict_key]['vds'],JFET[dict_key]['Id(J3)'],dict_key.split(' ')[2]+'V'])
            
            # generate resistance
            resistance = np.asarray(JFET[dict_key]['vds']) / np.asarray(JFET[dict_key]['Id(J3)'])
            
            # delete item 10
            resistance = np.delete(resistance,10)
            vds = np.delete(JFET[dict_key]['vds'],10)
            vgs = (dict_key.split(' ')[2]).replace('Vgs=','').replace('m','e-003')
            
            if vgs.find('-1'):
                if  vgs.find('-2'):
                    Plot1.append([vds, resistance,dict_key.split(' ')[2]+'V'])
                    
    for dict_key in JFET_lin.keys():
        
        # generate list
        if(len(dict_key.split(' ')) > 3):

            # generate list
            Plot_lin.append([JFET_lin[dict_key]['vds'],JFET_lin[dict_key]['Id(J3)'],(dict_key.split(' ')[2]+'V').replace('Vgs','Vatt')])
            
            # generate resistance
            resistance = np.asarray(JFET_lin[dict_key]['vds']) / np.asarray(JFET_lin[dict_key]['Id(J3)'])
            
            # delete item 10
            resistance = np.delete(resistance,10)
            vds = np.delete(JFET_lin[dict_key]['vds'],10)
            vgs = (dict_key.split(' ')[2]).replace('Vgs=','').replace('m','e-003')
            
            if vgs.find('-3'):
                if  vgs.find('-4'):
                    Plot1_lin.append([vds, resistance,(dict_key.split(' ')[2]+'V').replace('Vgs','Vatt')])   
            
            
#############################################################################
    # Plot FVCO
    Xlabel = ['$V_\mathrm{ds}$', 'V']
    Ylabel = ['$I_\mathrm{d}$', 'A']
       
    # Generate Plot
    plt.figure(figsize=(6,5))
    
    # some text
    ax1 = plt.subplot(111)
    pcb.Linear_Plot(ax1, Plot, Xlabel, Ylabel, Legend=True, LegendLoc = 4)
    ax1.set_ylim([-0.02,0.02])
    
    # Operation Region
    rect = patches.Rectangle((-1,-0.08),2,0.1,facecolor='#e6e6e6')
    ax1.add_patch(rect)
    ax1.annotate("", xy=(-1, 0.015), xytext=(1, 0.015), arrowprops=dict(arrowstyle="<->",linewidth=2))
    ax1.text(-0.05, 0.016, '$\overline{\mathrm{V}_\mathrm{\mathrm{A}-\overline{\mathrm{A}}}}$', fontsize=15)

    plt.savefig("JFET_Characteristic.pdf", dpi=120, bbox_inches='tight')
    
    # align axis
    plt.show   

    # Generate Plot
    plt.figure(figsize=(6,5))
    
    # some text
    ax1 = plt.subplot(111)
    pcb.Linear_Plot(ax1, Plot_lin, Xlabel, Ylabel, Legend=True, LegendLoc = 4)
    ax1.set_ylim([-0.02,0.02])
    
    # Operation Region
    rect = patches.Rectangle((-1,-0.08),2,0.1,facecolor='#e6e6e6')
    ax1.add_patch(rect)
    ax1.annotate("", xy=(-1, 0.015), xytext=(1, 0.015), arrowprops=dict(arrowstyle="<->",linewidth=2))
    ax1.text(-0.05, 0.016, '$\overline{\mathrm{V}_\mathrm{\mathrm{A}-\overline{\mathrm{A}}}}$', fontsize=15)

    # align axis
    plt.show     
    
    plt.savefig("JFET_Characteristic_lin.pdf", dpi=120, bbox_inches='tight')

#############################################################################
    # Plot FVCO
    Xlabel1 = ['$V_\mathrm{ds}$', 'V']
    Ylabel1 = ['$r_\mathrm{ds}$', '$\Omega$']
       
    # Generate Plot
    plt.figure(figsize=(6,5))
    
    # some text
    ax1 = plt.subplot(111)
    pcb.Linear_Plot(ax1, Plot1, Xlabel1, Ylabel1, Legend=True, LegendLoc = 2, marker='x')
    ax1.set_ylim([0,250])
    
        # Operation Region
    rect = patches.Rectangle((-1,-0.08),2,300,facecolor='#e6e6e6')
    ax1.add_patch(rect)
    
    # align axis
    plt.show       
    
    plt.savefig("JFET_Characteristic_RDS.pdf", dpi=120, bbox_inches='tight')            
           
 #############################################################################
    # Plot FVCO
    Xlabel1 = ['$V_\mathrm{ds}$', 'V']
    Ylabel1 = ['$r_\mathrm{ds}$', '$\Omega$']
       
    # Generate Plot
    plt.figure(figsize=(6,5))
    
    # some text
    ax1 = plt.subplot(111)
    pcb.Linear_Plot(ax1, Plot1_lin, Xlabel1, Ylabel1, Legend=True, LegendLoc = 2, marker='x')
    ax1.set_ylim([0,250])
    
        # Operation Region
    rect = patches.Rectangle((-1,-0.08),2,300,facecolor='#e6e6e6')
    ax1.add_patch(rect)
    
    # align axis
    plt.show       
    
    plt.savefig("JFET_Characteristic_RDS_lin.pdf", dpi=120, bbox_inches='tight')               

#############################################################################
    
if False:
    # Plot VCO Transfer Function
    XOR_SE = pcb.CSV2Dict("_JFET_Adder_Time.txt", delimiter="\t", blockheader=False, dictkey=True)
    Plot = []
    
    for dict_key in XOR_SE.keys():
        # generate list
        if(len(dict_key.split(' ')) > 3):
            # generate list
            label= (dict_key.split(' ')[2]+'V').replace('Vgs','Vatt')
            
            time = XOR_SE[dict_key]["time"][0:int(np.size(XOR_SE[dict_key]["time"])/3)]
            
            Plot.append([time,XOR_SE[dict_key]['V(out)'],label])

#############################################################################
    # Plot FVCO
    
    X_label = ['', 's']
    Ylabel_1 = [r'$\mathrm{V}_\mathrm{out}$', 'V']

    # Generate Plot
    plt.figure(figsize=(6,5))
    
    ax1 = plt.subplot(111)
    pcb.Linear_Plot(ax1, Plot, X_label, Ylabel_1)
    
    
    plt.show
    
    plt.savefig("JFET_SimpleCircuitry_TimeDomain.pdf", dpi=120, bbox_inches='tight')
   
 
    #############################################################################
    
if False:
    # Plot VCO Transfer Function
    XOR_SE = pcb.CSV2Dict("_JFET_Gain_transferfunction.txt", delimiter="\t", complexdelimiter='c')
    
#############################################################################
    # Plot FVCO
    X_label = [r'$\mathrm{V}_\mathrm{att}$', 'V']  
    Vatt = np.asarray(XOR_SE["v_att"])
    
    Resistance= np.asarray(XOR_SE["V(In,N001)/Id(J1)"])
    Linearized = Vatt*0.27+0.85
    
    [a,b,c] = pcb.FitFct_Exp(Vatt,Resistance, UsePoints='8:-1') 
    Fitted_Res = a + b*np.exp(c*Vatt)
    
    Ylabel_1 = [r'$\left| \frac{\mathrm{V}_\mathrm{out}}{\mathrm{V}_\mathrm{in}} \right|$', '']
    Plot_1 = [[Vatt, XOR_SE["-V(out)/V(in)"], "Gain"]]
    

    Plot_1_lin = [[Vatt,Linearized, "Linearized Gain"]]

    
    Ylabel_2 = [r'$\mathrm{r}_\mathrm{ds}$', '$\Omega$']
    Plot_2 = [[Vatt, Resistance, "$r_\mathrm{ds}$"]]
    Plot_2_Fit = [[Vatt, Fitted_Res, "approx. $r_\mathrm{ds}$"]]
    
    # Generate Plot
    plt.figure(figsize=(12,5))
    
    ax1 = plt.subplot(111)
    pcb.Linear_Plot(ax1, Plot_1, X_label, Ylabel_1)
    pcb.Linear_Plot(ax1, Plot_1_lin, X_label, Ylabel_1, Ylim=[0,1], linestyle='--')
    
    ax2 = ax1.twinx() 
    pcb.Linear_Plot(ax2, Plot_2, X_label, Ylabel_2, TwinX=ax1,
                    Ylim=[0,350], color='C2', LegendLoc=2)

    pcb.Linear_Plot(ax2, Plot_2_Fit, X_label, Ylabel_2, TwinX=ax1,
                    Ylim=[0,350], color='C3', LegendLoc=2, linestyle='--')    
    
    plt.show
    
    plt.savefig("JFET_SimpleCircuitry_TransferFunction.pdf", dpi=120, bbox_inches='tight')
    

 
    
#############################################################################
    
if True:
    # Plot VCO Transfer Function
    XOR_SE = pcb.CSV2Dict("_JFET_Gain_transferfunction.txt", delimiter="\t", complexdelimiter='c')
    XOR_SE_V12 = pcb.CSV2Dict("_JFET_Gain_transferfunction_V12.txt", delimiter="\t", complexdelimiter='c')   
    # measured parameters
    V_att = [0, -0.5, -1,-1.5, -2, -2.5, -3, -3.5, -4]
    V_out_A = [1.19, 0.954, 0.751, 0.572, 0.412, 0.32, 0.267, 0.251, 0.24]
    V_out_B = [1.12, 0.933, 0.719, 0.558, 0.4, 0.342, 0.282, 0.262, 0.252]
    
    Vin = 0.8
 
    V_out_A_Offset = [0.834, 0.628, 0.462, 0.312, 0.189, 0.102, 0.044, 0.018, 0.017]
    V_out_B_Offset = [0.819, 0.609, 0.430, 0.289, 0.164, 0.091, 0.037, 0.018, 0.013]
    
    G_A = np.multiply(V_out_A,1/Vin)
    G_B = np.multiply(V_out_B,1/Vin)

    G_A_Offset = np.multiply(V_out_A_Offset,1/Vin)
    G_B_Offset = np.multiply(V_out_B_Offset,1/Vin)
    
#############################################################################
    # Plot FVCO
    X_label = [r'$\mathrm{V}_\mathrm{att}$', 'V'] 
    Ylabel_1 = [r'$\left| \frac{\mathrm{V}_\mathrm{out}}{\mathrm{V}_\mathrm{in}} \right|$', '']
    
    Vatt = np.asarray(XOR_SE["v_att"])
    Vatt_12 = np.asarray(XOR_SE_V12["v_att"])
     
    Linearized = Vatt*0.27+0.85

    # Simulated
    Plot_1 = [[Vatt, XOR_SE["-V(out)/V(in)"], "sim. Gain (V1.0)"]]
    Plot_2 = [[Vatt_12, XOR_SE_V12["-V(out)/V(in)"], "sim. Gain (V1.2)"]]    
    
    # Approximated
    #Plot_1_lin = [[Vatt,Linearized, "sim. linearized Gain"]]
    
    # Measured
    #Plot_1_Meas =  [[V_att, G_A, "meas. Gain A (offset = 0V)"], [V_att, G_B, "meas. Gain B (offset = 0V)"],
     #                [V_att, G_A_Offset, "meas. Gain A (offset = -150mV)"], [V_att, G_B_Offset, "meas. Gain B (offset = -150mV)"]]
    
    # Generate Plot
    plt.figure(figsize=(12,5))
    
    ax1 = plt.subplot(111)
    pcb.Linear_Plot(ax1, Plot_1, X_label, Ylabel_1, linewidth=3)
    pcb.Linear_Plot(ax1, Plot_2, X_label, Ylabel_1, linewidth=3)
    #pcb.Linear_Plot(ax1, Plot_1_lin, X_label, Ylabel_1, linewidth=1.5, linestyle='--')
    #pcb.Linear_Plot(ax1, Plot_1_Meas, X_label, Ylabel_1, linewidth=1.5,  marker='x')    
  
    
    plt.show
    
    plt.savefig("JFET_SimpleCircuitry_TransferFunction_offset.pdf", dpi=120, bbox_inches='tight')
    

    
 #############################################################################   
    
if False:
    
    X_label = ['', 's']  
    # Plot VCO Transfer Function
    Adder_Complete = pcb.CSV2Dict("_ECL_LoopFilter_variable_K1.txt", delimiter="\t", complexdelimiter='c')
    
    # limit time
    time = Adder_Complete["time"][0:int(np.size(Adder_Complete["time"])*0.55)]

    #first plot
    Ylabel_1 = [r'$\overline{\mathrm{V}_\mathrm{\mathrm{x}-\overline{\mathrm{x}}}}$', 'V']
    Plot_1 = [[time, Adder_Complete["-V(a_se)"], "Input A"],
               [time, Adder_Complete["-V(b_se)"], "Input B"]]
    
    #second plot A
    Ylabel_2 = [r'$V_\mathrm{att,x}$', 'V']
    Plot_2 = [[time, Adder_Complete["V(att_a)"], "Attenuation Voltage A"],
               [time, Adder_Complete["V(att_b)"], "Attenuation Voltage B"]]
    
    #second plot B
    Atten_A = 0.27*np.asarray(Adder_Complete["V(att_a)"])+0.85
    Atten_B = 0.27*np.asarray(Adder_Complete["V(att_b)"])+0.85
   
    Ylabel_2A = [r'$G_k^\mathrm{a,1x}$', '']
    Plot_2A = [[time, Atten_A, "Gain A"],
               [time, Atten_B, "Gain B"]]
    
    # Calculated Channels
    Gain_A = Atten_A*np.asarray(Adder_Complete["-V(a_se)"])
    Gain_B = Atten_B*np.asarray(Adder_Complete["-V(b_se)"])
    Gain_AB = Gain_A+Gain_B
 
    Ylabel_3A = [r'$G_k^\mathrm{a,1x} \cdot \overline{\mathrm{V}_\mathrm{\mathrm{x}-\overline{\mathrm{x}}}}$', 'V']
    Plot_3A = [[time, Gain_A, "weighted input A"],
               [time, Gain_B, "weighted input B"]]

    Plot_3B = [[time, Gain_AB, "calculated Output"],
               [time, Adder_Complete["V(ab)"], "simulated Output"]]
    
    #third plot
    Ylabel_3 = [r'$V_\mathrm{add}$', 'V']
    Plot_3 = [[time, Adder_Complete["-V(b_se)-V(a_se)"], "A + B (w/o weighting)"],
               [time, Adder_Complete["V(ab)"], "A + B (w weighting)"]]
    
    #Vertical lines
    x_lines = [235e-9, 400e-9, 633e-9, 834e-9, 1035e-9]
    
    # Generate Plot
    plt.figure(figsize=(12,15))
    
    ax1 = plt.subplot(511)
    pcb.Linear_Plot(ax1, Plot_1, X_label, Ylabel_1, LegendLoc=2)

    ax2 = plt.subplot(512)
    pcb.Linear_Plot(ax2, Plot_2, X_label, Ylabel_2, LegendLoc=2)
    #pcb.Linear_Plot(ax2, Plot_2A, X_label, label_2A, Ylim=[-0.25,1.25], linestyle='--')

    ax2A = plt.subplot(513)
    pcb.Linear_Plot(ax2A, Plot_2A, X_label, Ylabel_2A, LegendLoc=2)
 
    ax3A = plt.subplot(514)
    pcb.Linear_Plot(ax3A, Plot_3A, X_label, Ylabel_3A, LegendLoc=2)
    
    ax3 = plt.subplot(515)
    pcb.Linear_Plot(ax3, Plot_3B, X_label, Ylabel_3, LegendLoc=2)
    
    
    # draw lines over everything
    for xline in x_lines:
        ax1.axvline(xline,ymin=-1,ymax=1.2 ,c="k",linewidth=2, clip_on=False, linestyle='--')
        ax2.axvline(xline,ymin=-1,ymax=1 ,c="k",linewidth=2, clip_on=False, linestyle='--')
        ax2A.axvline(xline,ymin=-1,ymax=1 ,c="k",linewidth=2, clip_on=False, linestyle='--')
        ax3A.axvline(xline,ymin=-1,ymax=1 ,c="k",linewidth=2, clip_on=False, linestyle='--')
        ax3.axvline(xline,ymin=0,ymax=1 ,c="k",linewidth=2, clip_on=False, linestyle='--')
        
    ax1.text(80e-9, 1.4, 'Section I', fontsize=12)
    ax1.text(280e-9, 1.4, 'Section II', fontsize=12)
    ax1.text(470e-9, 1.4, 'Section III', fontsize=12)  
    ax1.text(680e-9, 1.4, 'Section IV', fontsize=12)
    ax1.text(900e-9, 1.4, 'Section V', fontsize=12)        
    ax1.text(1050e-9, 1.4, 'Section VI', fontsize=12)
        
   
    
    plt.show
    
    plt.savefig("Adder_TransferExample.pdf", dpi=120, bbox_inches='tight')
    
#############################################################################
    
if False:
    # Plot VCO Transfer Function
    JFET_AC = pcb.CSV2Dict("_JFET_Gain_transferfunction_AC.txt", delimiter="\t", complexdelimiter=",",
                          blockheader=False, dictkey=True)
    JFET_AC.pop('')  
    
    # 3dB BW
    Plot = []
    
    # rename dict keys
    for dict_key in JFET_AC.keys():   
        # generate list
        if(len(dict_key.split(' ')) > 3):
            
            magnitue= np.abs(JFET_AC[dict_key]['V(out)/V(in)'])
            name = (dict_key.split(' ')[2]+'V').replace('V_att','$\mathrm{V}_\mathrm{att,x}$')
            # generate list
            Plot.append([JFET_AC[dict_key]['Freq.'],magnitue,name])
            
    JFET_AC_pos = pcb.CSV2Dict("_JFET_Gain_transferfunction_AC_posVGS.txt", delimiter="\t", complexdelimiter=",")
    Plot_pos = [[JFET_AC_pos['Freq.'],np.abs(JFET_AC_pos['V(out)/V(in)']),
                r'$\mathrm{V}_\mathrm{att,x}=\mathrm{1V}$']]
#############################################################################
    # Plot FVCO
    
    Xlabel = ['', 'Hz']
    Ylabel = [r'$G_k^{a,1x}$', '']
  
    # Generate Plot
    plt.figure(figsize=(12,5))
    
    ax1 = plt.subplot(111)
    pcb.SemiLogX_Plot(ax1, Plot, Xlabel, Ylabel, LegendLoc=4) 
    pcb.SemiLogX_Plot(ax1, Plot_pos, Xlabel, Ylabel, LegendLoc=4, linestyle='--') 
     
    # Operation Region
    rect = patches.Rectangle((1e6,-0.5),299e6,3,facecolor='#e6e6e6')
    ax1.add_patch(rect)
    
    ax1.annotate("", xy=(1e6, 1.2), xytext=(300e6, 1.2), arrowprops=dict(arrowstyle="<-",linewidth=2))
    ax1.text(12e6, 1.25, '$\mathrm{BW}_\mathrm{tune}$', fontsize=15)
    
    # align axis
    plt.show
    
    plt.savefig("JFET_Gain_AC.pdf", dpi=120, bbox_inches='tight')
    
#############################################################################
    
if False:
    # Plot VCO Transfer Function
    JFET_AC = pcb.CSV2Dict("HMC739_ModBW.txt", delimiter="\t", complexdelimiter=",",
                          blockheader=False, dictkey=True)
    JFET_AC.pop('')  
    Plot = []
    BW_Values = []
     
    BW_3dB_Y = -3
    
    for dict_key in JFET_AC.keys():   
        # generate list
        if(len(dict_key.split(' ')) > 3):
        
            # Plot
            dB_Plot = 10*np.log(np.abs(JFET_AC[dict_key]['V(n002)']))
            Freq_plot = JFET_AC[dict_key]['Freq.']
            name = (dict_key.split(' ')[2]+'F').replace('Vrac','$\mathrm{C}_\mathrm{var}$')
            
            Plot.append([Freq_plot,dB_Plot,name])
            
            # find 3dB Point
            BW_3dB = pcb.FindPoint_NextValue(dB_Plot, Freq_plot, BW_3dB_Y)
            BW_Values.append([BW_3dB,BW_3dB_Y,''])

        
#############################################################################
    # Plot FVCO  
    Xlabel = ['', 'Hz']
    Ylabel = [r'$G_\mathrm{VCO}$', 'dB']
  
    # Generate Plot
    plt.figure(figsize=(5,5))
    
    ax1 = plt.subplot(111)
    pcb.SemiLogX_Plot(ax1, Plot, Xlabel, Ylabel, Ylim=[-15,1], LegendLoc=3) 
    
    ax1.semilogx([10e6, 1e9], [-3,-3], 'k--')    
    pcb.SemiLogX_Plot(ax1, BW_Values, Xlabel, Ylabel,
                      XAutolim=False, LegendLoc=3, marker='H', color='r') 
    
    ax1.text(200e6, -2.5, "3dB-Bandwidth", bbox=dict(color='white', alpha=0.75))
     
    # Operation Region
    rect = patches.Rectangle((10e6,-20),289e6,25,facecolor='#e6e6e6')
    ax1.add_patch(rect)
    
    ax1.annotate("", xy=(10e6, -8), xytext=(300e6, -8), arrowprops=dict(arrowstyle="<-",linewidth=2))
    ax1.text(12e6, -7.5, '$\mathrm{BW}_\mathrm{tune}$', fontsize=15)
    
    # align axis
    plt.show
    
    plt.savefig("BW_VCOTune.pdf", dpi=120, bbox_inches='tight')    