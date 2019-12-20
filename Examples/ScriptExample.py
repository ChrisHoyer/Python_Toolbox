# Measure the PA output power - input power characterisitc using Aliglent_E8257D as input, power metre at the output.
# This script is modified by Mengqi from Zoltan's script. 
# Author: mengqi.cui@tu-dresden.de zoltan.tibenszky@tu-dresden.de
# Version: 08.07.2019
#
# Equipment: E8257D, 2ch HP6626A, +HP34401A DMM for current measurement


import visa, time, datetime, csv, logging, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import instruments.Instrument
import pandas as pd
from quantiphy import Quantity
from instruments.Agilent_E8257D import Agilent_E8257D
from instruments.RS_FSW import RS_FSW
from instruments.Agilent_6626A import Agilent_6626A
from instruments.Agilent_34410A import Agilent_34410A
from instruments.RS_ZVA67 import RS_ZVA67
from instruments.RS_NRP_Z55 import RS_NRP_Z55
###############################################################
def init_DC_source(dc_source, data_v, data_i):
    # dc_source.reset()
    for i in range(1,5):
        dc_source.set_vrange(voltage=1, channel=i)
        dc_source.set_current(current_limit=data_i[i-1], channel=i)
        dc_source.set_voltage(data_v[i-1], channel=i)
        dc_source.set_output_state(state=1, channel=i)

################################################################
def init_sg(sg):
    print("Initializing the Signal Generator")
    sg.set_output_disabled()
    sg.set_freq(freq)
    sg.set_output_power(-40) # 0dBm
    sg.set_output_enabled()

################################################################
def init_pm(pm):
    print("Initializing the Power Meter")
    pm.reset()
    #pm.set_zero()
    pm.set_freq(6.125e10)
    
####################################################################################################
def power_comp(freq, table):
    table_f = np.asarray([i[0] for i in table], dtype=np.float32)
    ind = np.abs(table_f - freq).argmin()
    return table[ind][1]

def dBm_to_W(x):
    return 10**(x/10)*1e-3
    
def pae(pin, pout, idd, vdd):
    return (dBm_to_W(pout)-dBm_to_W(pin))/(vdd*idd)

def drain_eff(pout, idd, vdd):
    return dBm_to_W(pout)/(idd*vdd)

####################################################################################################
# returns with a list of x,y pairs of the peak coordinates or {get_peak_list()} output ?FALSE? if noting
def meas_core_pm(fin, pin):
    # sg = signal generator instrument
    # fin: input freq
    # pin: input power
    # sa = spectrum analyser
    logging.info("Measuring freq=%s\tPin=%s" % (Quantity(fin,'Hz'),Quantity(pin,'dBm')))
    sg.set_output_power(pin) # - float(power_comp(freq, table=att_freq_data)))
    sg.set_freq(fin)
    pm.set_freq(fin)
    
    time.sleep(0.25) # just to be sure the SG changed the power

    power=pm.get_power()

    return [fin,power]

###
# Note: Bias Point:
# 1. VG1_0.25_0.525_VG2_1.1_VB2_VDD1.8
# 2. VG1_0.2_0.55_VG2_1.1_VB1.5_VDD1.8
# 3. VG1_0.3_0.55_VG2_1.1_VB1_VDD1.8
# 4. VG1_0.2_0.56_VG2_1.1_VB1.25_VDD1.8

        
if __name__ == '__main__':
    #######################################
    #   G L O B A L   V A R I A B L E S   #
    #######################################
    circuit_name   = "PA_V5_LS_55G_67G_Chip_A3_VG1_0.2_0.56_VG2_1.1_VB_1.25_VDD_1.8" 
    chip_number    = "1"
    freq           = 55e9
    fin_min        = 55e9
    fin_max        = 68e9
    fin_step       = 1e9
    pin_max        = 13   # dBm 
    pin_min        = -25  # dBm 
    pin_step       = 1  # dBm
    VDD         = 1.8
    plotting       = False
    measure        = False
    Rbt            = 3.07 
    
    [VG1_1,Nothing1 , VG1_2, Vdd] = [1,2,3,4] # DC source Betty GPIB27
    l_v_init_Betty = [0.2, VDD, 0.55, VDD]
    l_ilimit_Betty = ["1e-3", "500e-3", "1e-3", "500e-3"]
    
    [VG2, AM, VB, Nothing2] = [1,2,3,4] # DC source Barney   GPIB17
    l_v_init_Barney  = [1.1, 0, 1.25, 0]
    l_ilimit_Barney  = ["1e-3", "1e-3", "1000e-3", "1e-3"] 
    #dmm_Icore_Rmeas = 1 # 1Ohm

    logging.basicConfig(filename="PAV5_log_LS.log",level=logging.INFO)    
    logging.getLogger().addHandler(logging.StreamHandler())
    
    instruments.Instrument.logger.setLevel(logging.ERROR)
    visa.logger.setLevel(logging.ERROR)
    
    global fcsv_array
    fcsv_array = []
    
    # Instrument definition
    #dmm   = Agilent_34410A(addr="GBIP::22::INSTR")
    #sg    = Agilent_E8257D(addr= "TCPIP::172.31.228.97::INSTR")
    #sa    = RS_FSW(addr= "TCPIP::172.31.228.85::INSTR") # FSW
    Betty = Agilent_6626A(addr="GPIB1::27::INSTR",current_limit=1e-4, voltage_limit=0.9) # Betty
    Barney  = Agilent_6626A(addr="GPIB1::17::INSTR",current_limit=1e-4, voltage_limit=0.9) # Barney
    pm     = RS_NRP_Z55(addr="GPIB1::13::INSTR")
    sg     = Agilent_E8257D(addr= "GPIB1::20::INSTR")
    #ZVA = RS_ZVA67("TCPIP0::172.31.228.74::inst0::INSTR")
    #dmm   = Agilent_34410A(addr="TCPIP::172.31.228.115::INSTR")
    
    t_script_start = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    meas_date = time.strftime("%Y_%m_%d")   

    # Output file setup
    result_file = "%s_%s.csv" % (circuit_name, meas_date)
    log_file = result_file.replace(".csv",".log")


############################# Init ###################
    
    init_DC_source(Betty, l_v_init_Betty, l_ilimit_Betty)
    init_DC_source(Barney, l_v_init_Barney, l_ilimit_Barney)



    I_vdd=Betty.get_current(4)
    print(I_vdd)
    VDD_new=Rbt*I_vdd+VDD
    Betty.set_voltage(VDD_new,4)






    
    if measure:
        
        
        fp_csv = open(result_file,"w", newline='') # sensitivity curve data
        #fp_csv.write("freq[Hz],Pin_sg[dBm],Pin_fc[dBm],Pout_sa[dBm],Pout_fc[dBm],vdd[V],idd[A],PAE[],PE[]\n")
        fp_csv.write("freq[Hz],Pin_sg[dBm],Pout_Meas[dBm],I_DC[A],VDD_set[V],VDD\n")
        fcsv = csv.writer(fp_csv)

        init_sg(sg)
        init_pm(pm)

        fin=freq

        for fin in np.arange(fin_min,fin_max,fin_step):

            for pin in np.arange(pin_min,pin_max,pin_step):
                [fout, pout] = meas_core_pm(fin, pin)
                I_vdd=Betty.get_current(4)
                VDD_new=Rbt*I_vdd+VDD
                Betty.set_voltage(VDD_new,4)

                data = [fin, pin, pout,I_vdd, VDD_new, VDD]
                logging.debug(data)
                fcsv.writerow(data)

                
        
        # closing
        t_script_start = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        meas_date = time.strftime("%Y-%m-%d")
        meas_time = time.strftime("%H:%M")
        logging.info("Measurement finished at %s %s\n" % (meas_date, meas_time))

        Betty.set_voltage(0,4)
        Barney.set_voltage(0,1)
        Barney.set_voltage(0,3)
        sg.set_output_disabled()
        fp_csv.close()






    if plotting:
        # Plot 
        #   * Pout-Pin curve
        #   * VDD_set curve
        #   * ID curve
        df = pd.read_csv( result_file, header=0, index_col="freq[Hz]")
        print(df)
        pin  = df["Pin_sg[dBm]"].values
        pout = df["Pout_Meas[dBm]"].values
        idc = df["I_DC[A]"].values
        freq_meas=list(dict.fromkeys(df.index))
 
        
        fig,ax1 = plt.subplots()
        ax1.set_xlabel("Pin [dBm]")
        ax1.set_ylabel("Pout [dBm]")
        ax1.tick_params(axis="y")
        #ax1.set_ylim(-20,25)

        pinx=[]
        poutx=[]

        for i in range(len(freq_meas)):
            for j in range(len(df["Pin_sg[dBm]"])):
                if freq_meas[i] == df["Pin_sg[dBm]"].index[j]:
                    pinx.append(df["Pin_sg[dBm]"].values[j])
                    poutx.append(df["Pout_Meas[dBm]"].values[j])
            freqx=freq_meas[i]*1e-9
            ax1.plot(pinx, poutx)
            ax1.scatter(pinx, poutx, marker='s', s=10, label="%dGHz" % freqx)
            pinx=[]
            poutx=[]

            

            

        
        
        #ax2.set_xlabel("Pin [dBm]")
        #ax2.set_ylabel("I_DC", color=c2)
        #ax2.set_ylim(0,1)
        #ax2.plot(pin, idc,color=c2)
        #ax2.scatter(pin, idc,marker='s',s=10,label="I_DC",color=c2)
        #ax2.tick_params(axis="y",labelcolor=c2)
        
        #ax3 = ax1.twinx()
        #ax3.spines['right'].set_position(('outward', 60))
        #ax3.set_ylabel("PE [%]", color=c3)
        #ax3.plot(pin, pe,color=c3)
        #ax3.scatter(pin, pe,marker='o',s=10,label="pd",color=c3)
        #ax2.scatter(pin, pe,marker='o',s=10,label="pd",color=c3)
        fig.legend(loc="upper left", fancybox=True)
        ax1.grid()
        fig.tight_layout()
        #fig.savefig(result_file.replace(".csv",".png"), dpi=fig.dpi,bbox_inches='tight')
        plt.show()

