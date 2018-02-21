# =============================================================================
#  when frequencies match, outputs a uniform (ish) graph
# visually it is obvious when carrier and VCO frequency match, but must find a way to automate decision
# this copy amplifies the signal in each stage to make sure it is still maximum at 1.
# all variable names match diagram C.4.3 in GDP folder
# =============================================================================
import math as m
import scipy.signal as ss
import matplotlib.pyplot as plt
import numpy as np

VCO_START = 10000
SAMPLE_RATE = 480000
ENCODED_MESSAGE = "msk_encoded_message.txt"
MAX_SAMPLES = 10000
BIT_RATE = 6000.
BIT_PERIOD = 1/BIT_RATE


# =============================================================================
# Main functions
# =============================================================================

def frequency_finder(VCO_start=VCO_START, encoded_message=ENCODED_MESSAGE, sample_rate=SAMPLE_RATE):
    """
    Define parameters to find frequency.
    Calls a function that runs an iteration thrrough a costas loop
    Increments or decrements starting VCO_f depending on costas loop output
    Once conversion criteria are met, outputs converged VCO_f
    
    
    Variables:
        VCO_start               -   starting VCO frequency
        encoded_message         -   encoded message file input
        sample_rate             -   data sample rate
        
        encoded_message_file    -   file object to store encoded message file
        encoded_list            -   encoded message file amplitudes read and stored in list form
        f_found                 -   boolean to show if the carrier frequency has been found yet
        VCO_f                   -   variable to store current VCO frequency being tested
        
        costas_out              -   value quantifying if conversion has been met    
    """
    #open encoded message file and read samples to sample list converting to float
    encoded_message_file = open(encoded_message, "r")
    encoded_list = encoded_message_file.read().split("\n")
    
    #set butterworth filter variables
    B_cutoffAB2 = 0.03  
    B_orderAB2 = 5
    B_cutoffC3 = 0.01
    B_orderC3 = 5
    
    
    # three filter parameters,
    ABb, ABa = ss.butter(B_orderAB2, B_cutoffAB2, 'low')
    Cb, Ca = ss.butter(B_orderC3, B_cutoffC3, 'low')
    
    #conversion conditions
    f_found = False
    VCO_f = VCO_start
    
    # while conversion conditions are not met increment or decrement the VCO frequency
    while f_found == False:
        #costas out being the C channel term from costas loop diagram
        #VCO_f the current VCO frequeny being tested
        costas_out, VCO_f = costas_loop(VCO_f, encoded_list, ABa, ABb, Ca, Cb, sample_rate)
#        if costas_out == 0:
#            f_found = True
#        elif costas_out > 0:
#            VCO_f = VCO_f+1
#        elif costas_out < 0:
#            VCO_f = VCO_f-1
        f_found = True
    return VCO_f


def costas_loop(VCO_f, samples, filt1a, filt1b, filt2a, filt2b, sample_rate, max_samples=MAX_SAMPLES, T=BIT_PERIOD):
    #test plot variables
    frvC_A1_file = open("frvC_A1.txt", "w")
    frvC_B1_file = open("frvC_B1.txt", "w")
    frvC_C1_file = open("frvC_C1.txt", "w")
    frvC_A2_file = open("frvC_A2.txt", "w")
    frvC_B2_file = open("frvC_B2.txt", "w")
    frvC_A3_file = open("frvC_A3.txt", "w")
    frvC_B3_file = open("frvC_B3.txt", "w")
    frvC_C2_file = open("frvC_C2.txt", "w")
    frvC_C3_file = open("frvC_C3.txt", "w")
    frvC_C4_file = open("frvC_C4.txt", "w")
    timer_file = open("frvC_timer.txt", "w")
    
    #lists to store working variables
    ts = []
    A1 = []
    B1 = []
    C1 = []
    A2 = []
    B2 = []
    A3 = []
    B3 = []
    C2 = []
    C3 = []
    C4 = []
    
        
    #define timer variables
    t = 0
    sample_period = 1/sample_rate
    
    #oerform first multiplication: A1, B1
    for i in range(0, max_samples):
        ts.append(t)
        A1.append(float(samples[i])*m.sin(2*m.pi*VCO_f*t))
        B1.append(float(samples[i])*m.cos(2*m.pi*VCO_f*t))
        t += sample_period
        
        #generate timer file
        timer_file.write(str(t)+"\n")
    
    #pass through low pass filter: A2, B2
    A2 = ss.filtfilt(filt1b, filt1a, A1)
    B2 = ss.filtfilt(filt1b, filt1a, B1)
    
    # amplifiy A and B channels by 2
    for i in range(0, max_samples):
        A3.append(2*A2[i])
        B3.append(2*B2[i])
        
    #multiply A and B channels: C1
    for i in range(0, max_samples):
        C1.append(2*A3[i]*B3[i])

    #multiply by sin((pi/T)t: C2
    for i in range(0, max_samples):
        C2.append(C1[i]*m.cos((m.pi/T)*ts[i]))
    
    
    # pass through low pass filter: C3
    C3 = ss.filtfilt(filt2b, filt2a, C2)
     
    #square so all values are positive
    for i in range(0, max_samples):
        C4.append(C3[i]**2)
    
    #integrate over entire samples, if perfectly matched and perfectly filtered integration will equal zero
    f_err = sum(C4)
    
    print(f_err)
    #writing to test files
    for i in range(0, max_samples):
        frvC_A1_file.write(str(A1[i])+"\n")
        frvC_B1_file.write(str(B1[i])+"\n")
        frvC_A2_file.write(str(A2[i])+"\n")
        frvC_B2_file.write(str(B2[i])+"\n")
        frvC_A3_file.write(str(A3[i])+"\n")
        frvC_B3_file.write(str(B3[i])+"\n")
        frvC_C1_file.write(str(C1[i])+"\n")
        frvC_C2_file.write(str(C2[i])+"\n")
        frvC_C3_file.write(str(C3[i])+"\n")
        frvC_C4_file.write(str(C4[i])+"\n")
    
    #close all open files
    frvC_A1_file.close()
    frvC_B1_file.close()
    frvC_C1_file.close()
    frvC_A2_file.close()
    frvC_B2_file.close()
    frvC_A3_file.close()
    frvC_B3_file.close()
    frvC_C2_file.close()
    frvC_C3_file.close()
    frvC_C4_file.close()

    return VCO_f, f_err


# =============================================================================
# Other  functions for debugging
# =============================================================================
    
def biggest_magnitude(inlist):
    """
    Takes in a list and return value with largest magnitude regardless of sign
    """
    abslist = []
    for i in inlist:  
        if i < 0:
            abslist.append(i*-1)
        else:
            abslist.append(i)
    return inlist[abslist.index(max(abslist))]

def is_equal(inlist1, inlist2):
    """
    Compares two lists and outputs True if they are identical and False if they are not
    """
    equal = True
    discr_elements = []
    if len(inlist1) != len(inlist2):
        print("list lengths not equal")
        return False
    else:
        for i in range(0, len(inlist1)):
            if inlist1[i] != inlist2[i]:
                equal = False
                discr_elements.append(i)
    
    print(equal)
    print("Discrepency elements: ")
    print(discr_elements)
    return equal, discr_elements

def frequency_sweep(f_start, f_end, encoded_message=ENCODED_MESSAGE, sample_rate=SAMPLE_RATE):
    """
    sweeps through a range of frequencies and quanitfies the error for each then graphs err against frequency
    """
    #open encoded message file and read samples to sample list converting to float
    encoded_message_file = open(encoded_message, "r")
    encoded_list = encoded_message_file.read().split("\n")
    
    #set butterworth filter variables
    B_cutoffAB = 0.03  
    B_orderAB = 5
    B_cutoffC = 0.5
    B_orderC = 1
    
    # three filter parameters,
    ABb, ABa = ss.butter(B_orderAB, B_cutoffAB, 'low')
    Cb, Ca = ss.butter(B_orderC, B_cutoffC, 'low')
    
    # variables to store plot data    
    frequency_list = []
    err_list = []
    
    for f in range(f_start, f_end):
        VCO_f, VCO_err = costas_loop(f, encoded_list, ABa, ABb, Ca, Cb, sample_rate)
        frequency_list.append(f)
        err_list.append(VCO_err)
    
    # plot the errors against frequency
    plt.plot(frequency_list, err_list, label="error")
    plt.show()
    
    return VCO_f
            
# =============================================================================
# Running code 
# =============================================================================
frequency_finder()
#frequency_sweep(9000, 11000)