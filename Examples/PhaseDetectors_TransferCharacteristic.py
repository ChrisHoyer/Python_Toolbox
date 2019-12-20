# packages
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import numpy as np

############################################################################# 
# number of points 
num_points = 1e3
# freq constant
omega0 = 1

# Phase
phase = np.linspace(-4*np.pi,4*np.pi,num_points)
time = np.linspace(0,num_points,num_points)  

# Two Signals, one is phase shifted
sref = np.sqrt(2)*np.cos(omega0*time)                      #w = fixed frequency
sfb = np.sqrt(2)*np.cos(omega0*time+phase) 


# Plot Input Signals
if False:
    plt.figure(figsize=(10, 5))
    plt.plot(phase, sfb, label="FB")
    plt.plot(phase, sref, label="REF")
    plt.grid()
    plt.legend()
    plt.show

############################################################################# 
# Lowpass Filter for Signal (rect)
def AVG_TF(input_signal, cutoff=0.1):
    
    fft_signal = np.fft.fft(input_signal)
    
    # length
    signal_length = np.size(fft_signal)
    filter_length = int(signal_length*cutoff)
    
    # lowpass
    passband = np.ones(filter_length)
    stopband = np.zeros(signal_length-filter_length)
    signal_filter = np.append(passband, stopband) 
    
    # use filter
    filtered_signal = fft_signal*signal_filter
    
    output_signal = np.fft.ifft(filtered_signal)
    
    # Plot Filter
    if False:
        plt.figure()
        plt.plot(np.real(fft_signal), 'b')
        plt.plot((np.real(signal_filter)*np.max(fft_signal)), 'k--')       
        plt.plot(np.real(filtered_signal), 'r')
        plt.show
        
    return output_signal

# Schmitt Trigger
def SchmittTrigger(input_signal, threshold=0, high=1, low=-1, edge=False):
    
    # Map for digital signals
    output_signal = input_signal-np.mean(input_signal)
    output_signal = np.where(output_signal > threshold, high ,low)
    
    if edge:

        # shift to find edge
        output_signal_shift = np.roll(output_signal,1)
        
        #edge detection
        output_signal = output_signal-output_signal_shift
        
        # shift back and remove extended
        output_signal = np.roll(output_signal,-1)
        
        #only rising edge
        output_signal = np.where(output_signal < threshold, high ,low)
    

    # Plot Schmitt Trigger
    if False:
        plt.figure()
        plt.plot(input_signal, 'r--')
        plt.plot(output_signal,'b')
        plt.show
    
    return output_signal
 
    
#############################################################################

# Function PD Multiply, Analog
def PD_Multiply(ref, fb):
    # Simple Multiplication
    PD_out = ref*fb
    return AVG_TF(PD_out, cutoff=0.01)


# Function with XOR
def PD_XOR(ref, fb):
    # Schmitt Trigger Input
    ref_sqrt = SchmittTrigger(ref, high=1, low=0)
    fb_sqrt = SchmittTrigger(fb, high=1, low=0) 
    
    # Bitwise XOR  # reduce DC Offset and scale by 2
    PD_out = np.bitwise_xor(fb_sqrt, ref_sqrt)
    PD_out = 2*(PD_out - np.mean(PD_out))
           
    return AVG_TF(PD_out, cutoff=0.00)

# Function with XOR
def PD_XOR_IDEAL(phase):
    
    # create the changing points
    PD_out_x = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    PD_out_y = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5]
    
    PD_out  = [PD_out_x, PD_out_y]
    
    return PD_out

# Function with D-FlipFlop (ideal)
def PD_DFF(ref, fb):
    
    # Up and Down Values
    Up_out = []
    Dwn_out = []
        
    # Schmitt Trigger Input with rising edge detection
    A_sqrt_edge = SchmittTrigger(fb, edge=True, high=1, low=0)
    B_sqrt_edge = SchmittTrigger(ref, edge=True, high=1, low=0) 
    
    # starting state
    Up = 0
    Down = 0
    
    # Iterate all items
    for i in range(np.size(A_sqrt_edge)):
        
        # State Diagramm - Page 5
        #http://pallen.ece.gatech.edu/Academic/ECE_6440/Summer_2003/L070-DPLL(2UP).pdf
        
        # rising edge signal A and Latch in reset
        if (A_sqrt_edge[i] == 1) & (Up == 0):
            Up = 1 
                
        # rising edge signal B  and Latch in reset              
        if (B_sqrt_edge[i] == 1) & (Down == 0):
            Down = 1
        
        # Switch in the middle to finally get a result...
        if i == int(np.size(A_sqrt_edge)/2):
            Up = 1
            Down = 0
            print("\n edge at " + str(i) + "\n")
        
        # Both Latches are in set state -> reset
        if (Up == 1) & (Down == 1):
            Up = 0
            Down = 0
                       
        Up_out.append(Up)
        Dwn_out.append(Down)
    
    # construct output signal
    PD_out = -1*(np.asarray(Up_out) - np.asarray(Dwn_out))
    
    # Plot Input Signals
    if True:
        middle_span = 50
        middle = np.size(A_sqrt_edge)/2
        
        f = plt.figure(figsize=(15,10))
        ax1 = f.add_subplot(311)
        ax2 = f.add_subplot(312)
        ax3 = f.add_subplot(313)
        #ax1.plot(A_sqrt, 'IndianRed', label="A (Up)")
        #ax1.plot(B_sqrt, 'SkyBlue', label="B (Down)")
        ax1.stem(A_sqrt_edge, linefmt='IndianRed', markerfmt=" ")
        ax1.stem(B_sqrt_edge, linefmt='SkyBlue', markerfmt=" ")
        ax1.grid()
        ax1.set_xlim(middle-middle_span, middle+middle_span)
        ax1.legend()
        ax2.plot(Up_out, 'IndianRed', label="Up")
        ax2.plot(Dwn_out, 'SkyBlue', label="Down")
        ax2.grid()
        ax2.set_xlim(middle-middle_span, middle+middle_span)
        ax2.legend()
        ax3.plot(PD_out, 'Teal', label="PD")
        ax3.grid()
        ax3.legend()
        ax3.set_xlim(middle-middle_span, middle+middle_span)
        plt.show 
        
    #return PD_out
    return AVG_TF(PD_out, cutoff=0.02)

# Function with XOR
def PD_DFF_IDEAL(phase):
    
    # create the changing points
    PD_out_x = [-4, -1.99, -2.01, 1.99,2.01,4]
    PD_out_y = [-0.5, 0, -0.5, 0.5,0,0.5 ]
    
    PD_out  = [PD_out_x, PD_out_y]
    
    return PD_out

#############################################################################

# Calculate Transfer Fuction
PD_ANALOG = PD_Multiply(sref, sfb)
PD_XORGate = PD_XOR_IDEAL(phase)
PD_DLatch = PD_DFF_IDEAL(phase)

#############################################################################
phase_pi = phase/np.pi
  
  
f,ax = plt.subplots(figsize=(10, 4))
ax.plot(phase_pi, PD_ANALOG, label="Analog PD")
ax.plot(PD_XORGate[0],PD_XORGate[1], label="XOR PD")
plt.plot(PD_DLatch[0],PD_DLatch[1], label="D-FlipFlop PFD")

# Labeling x-axis
ax.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
ax.xaxis.set_major_locator(plt.MultipleLocator(1))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))

# Labeling y-axis
ax.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $V_{out}$'))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

# grid
plt.minorticks_on()
plt.grid(True, which='major')
ax.grid(which='minor', alpha=1, linestyle=':', linewidth=1)
ax.grid(which='major', alpha=1, linewidth=1.2)

# labeling
ax.set_xlabel('$\Delta~\Phi$')
ax.legend(loc=2)
plt.show

plt.savefig("PD_transferfunction.pdf", dpi=120, bbox_inches='tight')


