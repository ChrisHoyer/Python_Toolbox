# Measure PVCO output power vs Vtune and Pdc
# This script is modified by Mengqi from Zoltan's script. 
# Author: mengqi.cui@tu-dresden.de zoltan.tibenszky@tu-dresden.de
# Version: 07.10.2019
#
# Equipment: 2X HP6626A, RS_FSU67 the RS_FSW.py is here used for this FSU67


import visa, time, datetime, csv, logging, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import instruments.Instrument
import pandas as pd
from quantiphy import Quantity
from instruments.Agilent_6626A import Agilent_6626A
from instruments.RS_NRP_Z55 import RS_NRP_Z55
from instruments.RS_FSW import RS_FSW


###############################################################
def init_DC_source(dc_source, data_v, data_i):
    # dc_source.reset()
    for i in range(1,5):
        dc_source.set_vrange(voltage=1, channel=i)
        dc_source.set_current(current_limit=data_i[i-1], channel=i)
        dc_source.set_voltage(data_v[i-1], channel=i)
        dc_source.set_output_state(state=1, channel=i)

################################################################
def init_sa(sa):
    print("Initializing the Spectrum Analyser")
    sa.reset()
    #sa.set_single_sweep_mode()
    sa.set_freq_center(freq)
    sa.set_freq_span(4e9)
    sa.set_resolution_bw(10e6)       # res bandwidth has a influence on pout measurement
    sa.set_sweep_points(1001)
    sa.set_unit("DBM")
    sa.set_avg()
    sa.set_avg_cnt(10)
    sa.move_marker_peak(1)
    sa.set_marker_count()
    sa.set_marker_peak_excursion(10)
    #sa.start_peak_search()
    #sa.set_marker_threshold_state("ON")
    #sa.set_marker_threshold(-40)
    #sa.sort_peak_list("Y")
    sa.set_display_update()
####################################################################################################


def dBm_to_W(x):
    return 10**(x/10)*1e-3
def dBm_to_mw(x):
    return 10**(x/10)


####################################################################################################


####################################################################################################
# returns with a list of x,y pairs of the peak coordinates or {get_peak_list()} output ?FALSE? if noting
# modifiey for Pvco measurement
def meas_core_vco(vtune):

    #logging.info("Measuring vtune=%s\tPin=%s" % (Quantity(fin,'Hz'),Quantity(pin,'dBm')))
    #sg.set_output_power(pin) # - float(power_comp(freq, table=att_freq_data)))
    #sg.set_freq(fin)
    #sa.set_freq_center(fin)
    wilma.set_voltage(vtune,3)
    fsu.set_freq_start(50e9)
    fsu.set_freq_stop(67e9)
    time.sleep(2) # just to be sure the SG changed the power

    fsu.move_marker_peak(1)     # the X and Y coordinate will be updated after about 1 sec
    time.sleep(1.5)
    yyy=fsu.get_marker_y(1)
    xxx=fsu.get_marker_x(1)
    #print(xxx)
    #print(yyy)
    fsu.set_freq_center(xxx)
    fsu.set_freq_span(1500e6)
    
    time.sleep(1.5)
    fsu.move_marker_peak(1)
    time.sleep(1.5)
    yyy=fsu.get_marker_y(1)
    xxx=fsu.get_marker_x(1)
    #print(xxx)
    #print(yyy)
    fsu.set_freq_center(xxx)
    fsu.set_freq_span(500e6)

    fsu.single_run()
    fsu.start_peak_search() # sorted by power
    pl = fsu.get_peak_list()#
    #print(pl)
    if pl:
        logging.debug(pl)
        f_list = [i[0] for i in pl]
        p_list = [i[1] for i in pl]
        cnt = 0
        max = 0
        for i in range(len(p_list)):
            if p_list[i] > max:
                max = p_list[i]
                cnt = i
        return [f_list[cnt],p_list[cnt]]        

    else: # pl=False --> No peaks detected
        logging.debug("No peak found.")
        return False

####################################################################################################



    #########################################################
    #  Start !!                                             #
    #########################################################

       
if __name__ == '__main__':
    #######################################
    #   G L O B A L   V A R I A B L E S   #
    #######################################
    circuit_name   = "PVCO_ChipA2_Vara_VDDCore_1p_VGC_0p525_VDDBuf_2p7_VG12_0p475_1p35" 
    chip_number    = "1"
    freq           = 64e9
    #pin            = 0 #dBm   
    #fin_min        = 55e9
    #fin_max        = 67e9
    #fin_step       = 0.1e9
    #pin_max        = 20   # dBm 
    #pin_min        = -20  # dBm 
    #pin_min        = -20  # dBm
    pin_step       = 1  # dBm
    # pin_step       = 5  # dBm
    vdd_buffer     = 2.2
    vdd_core       = 1
    vg_core        =0.525
    vtune          = 1
    Rdc_line       = 0.5+0.25+0.2
    plotting       = False
    measure_vco     = True
    test           = False
    

    [ch_VDDbuf, ch_G2, ch_G1, ch_NA1] = [1,2,3,4] # DC source Betty
    l_v_init_betty  = [vdd_buffer, 1.35, 0.475, 0]
    l_ilimit_betty  = ["250e-3", "1e-3", "1e-3", "1e-3"]     
    [ch_CoreVG, ch_VDDCore, ch_vtune, ch_NA2] = [1,2,3,4] # DC source WILMA II
    l_v_init_wilma = [vg_core, vdd_core, vtune, 0]
    l_ilimit_wilma = ["1e-3", "200e-3", "1e-3", "1e-3"] 

  

    logging.basicConfig(filename="PVCO_Measurement.log",level=logging.INFO)    
    logging.getLogger().addHandler(logging.StreamHandler())
    
    instruments.Instrument.logger.setLevel(logging.ERROR)
    visa.logger.setLevel(logging.ERROR)
    
    global fcsv_array
    fcsv_array = []
    
    # Instrument definition
    # pm     = RS_NRP_Z55(addr="GPIB1::13::INSTR")
    fsu    = RS_FSW(addr= "GPIB1::20::INSTR") # FSU using the FSW.py, should insert more command if needed
    wilma = Agilent_6626A(addr="GPIB1::3::INSTR",current_limit=1e-4, voltage_limit=0.9) # Wilma II GPIB 3
    betty  = Agilent_6626A(addr="GPIB1::27::INSTR",current_limit=1e-4, voltage_limit=0.9) # Betty GPIB 27
   
    
    t_script_start = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    meas_date = time.strftime("%Y_%m_%d")
    

    # Output file setup
    result_file = "%s_%s.csv" % (circuit_name, meas_date)
    log_file = result_file.replace(".csv",".log")



    if test:
        #init_sg(sg)
        #init_pm(pm)
        init_sa(fsu)
        init_DC_source(wilma, l_v_init_wilma, l_ilimit_wilma)
        init_DC_source(betty, l_v_init_betty, l_ilimit_betty)
        
        res=meas_core_vco(1)
        print(res)




       # fsu.set_resolution_bw(1e5)

       # for i in range(5):
       #     time.sleep(0.1)
       #     x=5-i
       #     print("Shut down DC source in %d seconds." % x)

        #betty.set_output_state(0,"all")
        #wilma.set_output_state(0,"all")
            
        #meas_core_vco(vtune)
        #print("Measuring at Vtune = %d." % vtune)    
        









    
    if measure_vco:
        
        
        fp_csv = open(result_file,"w", newline='') # sensitivity curve data
        #fp_csv.write("freq[Hz],Pin_sg[dBm],Pin_fc[dBm],Pout_sa[dBm],Pout_fc[dBm],vdd[V],idd[A],PAE[],PE[]\n")
        fp_csv.write("Vtune[V],Freq[Hz],Pout_Meas[dBm],VDDbuff[V],Ibuff[mA],VDDCore[V],ICore[mA],Effi_raw\n")
        fcsv = csv.writer(fp_csv)


        #with open("sn06149358+mt77+sn10347319_slow.csv",newline='') as csvfile:
        #    att_freq_data = list(csv.reader(csvfile))[1:] 
            
        #####################################################
        #   E N D   O F   G L O B A L   V A R I A B L E S   #
        #####################################################

        init_DC_source(wilma, l_v_init_wilma, l_ilimit_wilma)
        init_DC_source(betty, l_v_init_betty, l_ilimit_betty)
        init_sa(fsu)

        Ibuff_re = betty.get_current(1) # A
        Icore_re = wilma.get_current(2)

        VDDbuff_new = vdd_buffer + Ibuff_re*Rdc_line
        VDDcore_new = vdd_core + Icore_re*Rdc_line

        betty.set_voltage(VDDbuff_new,1)
        wilma.set_voltage(VDDcore_new,2)

        

        vtune_list = np.arange(0, 2, 0.1)
        print("Measuring at VDDbuff: %0.2f, VDDCore: %0.2f, VG_core: %0.3f" %(vdd_buffer,vdd_core,vg_core))

        for vx in vtune_list:
            Ibuff_re = betty.get_current(1) # A
            Icore_re = wilma.get_current(2)
            VDDbuff_new = vdd_buffer + Ibuff_re*Rdc_line
            VDDcore_new = vdd_core + Icore_re*Rdc_line
            betty.set_voltage(VDDbuff_new,1)
            wilma.set_voltage(VDDcore_new,2)
            time.sleep(0.5) # wait for the voltage change
            
            res=meas_core_vco(vx)
            osc_freq = res[0]
            pout = res[1]
            Ibuff = betty.get_current(1)*1000 # mA
            Icore = wilma.get_current(2)*1000 
            Eff_raw = dBm_to_mw(pout)/(0.5*vdd_buffer*Ibuff+0.5*vdd_core*Icore)
            data = [vx,osc_freq,pout,vdd_buffer,vdd_core,Eff_raw]
            print("Measuring: vtune: %0.3f, osc_freq: %0.2f, Pout: %0.4f, Effi_raw: %0.4f" %(vx, osc_freq, pout, Eff_raw))
            logging.debug(data)
            fcsv.writerow(data)
        
        # closing
        t_script_start = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        meas_date = time.strftime("%Y-%m-%d")
        meas_time = time.strftime("%H:%M")
        logging.info("Measurement finished at %s %s\n" % (meas_date, meas_time))
        wilma.set_voltage(0,3)
        fp_csv.close()









    if plotting:
        # Plot 
        #   * Loss vs Freq curve
        #   
        #   
        #
        #result_file="AgilPSG_FSW_1-67G_Cable_SN11167431_fe2fe_2019_07_09.csv"
        
        df = pd.read_csv( result_file, header=0, index_col="freq[Hz]")
        print(df)
        
        pin  = df["Pin_sg[dBm]"].values
        pout = df["Pout_Meas[dBm]"].values

        line_loss = []
        for i in range(len(pin)):
            loss=pin[i]-pout[i]
            line_loss.append(loss)

        
        fig,ax1 = plt.subplots()
        c1 = "black"
        c2 = "red"
        c3 = "blue"
        
        ax1.set_xlabel("Freq [Hz]")
        ax1.set_ylabel("Line Loss [dB]", color=c1)
        ax1.set_ylim(-3,10)
        ax1.plot(df.index, line_loss, color=c1)
        ax1.scatter(df.index, line_loss, marker='s', s=10, label="Line_Loss",color=c1)
        ax1.tick_params(axis="y",labelcolor=c1)
        
        #ax2 = ax1.twinx()
        #ax2.set_ylabel("PAE & PE [%]", color=c2)
        #ax2.set_ylim(0,14)
        #ax2.plot(pin, pae,color=c2)
        #ax2.scatter(pin, pae,marker='+',s=10,label="pae",color=c2)
        #ax2.tick_params(axis="y",labelcolor=c2)
        
        #ax3 = ax1.twinx()
        #ax3.spines['right'].set_position(('outward', 60))
        #ax3.set_ylabel("PE [%]", color=c3)
        #ax3.plot(pin, pe,color=c3)
        #ax3.scatter(pin, pe,marker='o',s=10,label="pd",color=c3)
        #ax2.scatter(pin, pe,marker='o',s=10,label="pd",color=c3)
        #fig.legend(loc="upper left", fancybox=True)
        ax1.grid()
        fig.tight_layout()
        #fig.savefig(result_file.replace(".csv",".png"), dpi=fig.dpi,bbox_inches='tight')
        plt.show()

