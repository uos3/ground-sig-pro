# =============================================================================
# v1.6
# MSK working dependant on knowing carrier frequency
# - still contains writes to output files showing singal through processing
# - uses integration over symbol period to decided symbol
# - syncs to start of byte before outputting bits to output file
# ==============================

import math as m
import numpy as np
import scipy.signal as ss

CARRIER_FREQUENCY = 10000     # Hz
CARRIER_AMPLITUDE = 1       # V?
CARRIER_TIMESPAN = 5        # s
DC_SHIFT = 0                # V?
SAMPLE_RATE = 480000      # Hz
BIT_RATE = 6000              # bps
m_phases = 2
B_ORDER1 = 5
B_CUTOFF1 = 0.02
MESSAGE = "binary.txt"
ENCODED_MESSAGE = "msk_encoded_message.txt"
SYNC_KEY = '00101010'
B_ORDER = 5
B_CUTOFF = 0.05
    

def main():
    """
    Menu function to select what to do
    """

    print("1. msk encode\n2. msk decode")
    user_input = input("->")

    if user_input == "1":
        msk_encode()
    elif user_input == "2":
        msk_decode()
    elif user_input == "3":
        generate_carrier()
    elif user_input == "exit":
        pass
    else:        
        main()
    return 0

# Functions used in mainu
def msk_encode(carrier_A = CARRIER_AMPLITUDE, carrier_f = CARRIER_FREQUENCY, bit_rate=BIT_RATE, sample_rate=SAMPLE_RATE, message = MESSAGE):
    message_file = open(MESSAGE, "r")
    output_file = open("msk_encoded_message.txt", "w")
    mibits = open("mibits.txt", "w")
    mqbits = open("mqbits.txt", "w")
    mibits_msk = open("mibits_msk.txt", "w")
    mqbits_msk = open("mqbits_msk.txt", "w")
    miarm = open("miarm.txt.", "w")
    mqarm = open("mqarm.txt", "w")
    
    message = [i for i in message_file.read()]
    samples_per_bit = sample_rate/bit_rate
    samples_per_symbol = int(m_phases*sample_rate/bit_rate)
    sample_interval = 1/sample_rate
    bit_period = 1/bit_rate
    
    #convert to NRZ
    for i in range(0, len(message)): 
        message[i] = 2*int(message[i]) - 1
    
    #split odd and even bits into I bitstream and Q bitstream
    I_bits = message[0::2]
    Q_bits = message[1::2]
    
    #set variables for output
    t = -bit_period
    I_arm = []
    Q_arm = []
    sum_channel = []
    
    #for bit pair
    for i in range(0, len(I_bits)):
        #perform shift for as many samples in bit pair work out I and Q channels
        for sample in range(0, samples_per_symbol):
            
            # fixes fact that Q arm is on one bit period offset
            if t < 0:
                tQ = 0
            else:
                tQ = t
            
            # fixes changing what bit Q arm is on mid bitform
            if sample < samples_per_symbol/2:
                Qi = i-1
            else:
                Qi = i
                
            #for visualisation and testing of I and Q bitstreams
            mibits.write(str(I_bits[i])+"\n")
            mqbits.write(str(Q_bits[Qi])+"\n")
            
            #for visualisation an testing of MSK I and Q bits forms
            mibits_msk.write(str(I_bits[i]*m.cos((m.pi*t)/(2*bit_period)))+"\n")
            mqbits_msk.write(str(Q_bits[Qi]*m.sin((m.pi*tQ)/(2*bit_period)))+"\n")
            
            
            #perform I and Q arm modulation
            I_arm.append(I_bits[i]*m.cos((m.pi*t)/(2*bit_period))*m.cos(2.0*m.pi*carrier_f*t))
            Q_arm.append(Q_bits[Qi]*m.sin((m.pi*tQ)/(2*bit_period))*m.sin(2.0*m.pi*carrier_f*t))           
            
            #increment timer
            t += sample_interval
            
    # sum channels for final output value
    for i in range(0, len(I_arm)):
        sum_channel.append(I_arm[i] - Q_arm[i])
    
    #normalise output
    maxsc = max(sum_channel)
    for i in range(0, len(sum_channel)):
        sum_channel[i] = sum_channel[i]/maxsc
    
    #write output to output file
    for i in range(0, len(sum_channel)):
        output_file.write(str(sum_channel[i])+"\n")
        #for visualisation anf testing testing of I and Q modulated signals
        miarm.write(str(I_arm[i])+"\n")
        mqarm.write(str(Q_arm[i])+"\n")
    
    
    #close all files
    message_file.close()
    output_file.close()
    mibits.close()
    mqbits.close()
    mibits_msk.close()
    mqbits_msk.close()
    miarm.close()
    mqarm.close()
        
    return 0


def msk_decode(sample_rate = SAMPLE_RATE, bit_rate = BIT_RATE, encoded_filename = ENCODED_MESSAGE, sync_key=SYNC_KEY, B_order=B_ORDER, B_cutoff=B_CUTOFF):
#   open input and output files
    encoded_message = open("msk_encoded_message.txt", "r")
    output_file = open("msk_decoded_message.bin", "w")
    mdi1 = open("mdi1.txt.", "w")
    mdq1 = open("mdq1.txt", "w")
    mdi2 = open("mdi2.txt", "w")
    mdq2 = open("mdq2.txt", "w")
    mdi3 = open("mdi3.txt", "w")
    mdq3 = open("mdq3.txt", "w")
    mdi4 = open("mdi4.txt", "w")
    mdq4 = open("mdq4.txt", "w")
    mdibits = open("mditbits.txt", "w")
    mdqbits = open("mdqbits.txt", "w")
    
#   define constants
#   frequency recovery loop goes here
    carrier_f = CARRIER_FREQUENCY
    sample_period = 1/sample_rate
    samples_per_symbol = m_phases*sample_rate/bit_rate
    bit_period = 1/bit_rate
    
#   butterworth filter coefficients
    b, a = ss.butter(B_order, B_cutoff, 'low')
    
#   convert message to list
    encoded_list = encoded_message.read().split("\n")
    
    for index in range(0, len(encoded_list)):
        try:
            encoded_list[index] = float(encoded_list[index])
        except ValueError:
            del encoded_list[index]

    
    # define variables to hold output data
    I1 = []
    I2 = []
    I3 = []
    I4 = []
    I5 = []
    Q1 = []
    Q2 = []
    Q3 = []
    Q4 = []
    Q5 = []
    output = []
    sync_list = []
        
    t = 0
    message_length = len(encoded_list)
    message_sync = False
    
    for index in range(0, message_length):
        I1.append(encoded_list[index]*m.cos(2*m.pi*carrier_f*t))
        Q1.append(encoded_list[index]*m.sin(2*m.pi*carrier_f*t))
        
#       output for visualisation during testing
        mdi1.write(str(I1[index])+"\n")
        mdq1.write(str(Q1[index])+"\n")
        
        I2.append(I1[index]*m.cos((m.pi*t)/(2*bit_period)))
        Q2.append(Q1[index]*m.sin((m.pi*t)/(2*bit_period)))

#       output for visualisation during testing
        mdi2.write(str(I2[index])+"\n")
        mdq2.write(str(Q2[index])+"\n")
        
        t += sample_period

#   add shift to q list
    for time_period in np.arange(0, bit_period, sample_period):
        I2.append(0)
        Q2.insert(0, 0)
    
    message_length = len(I2)
    
#   for visualisation during testing
    for i in range(0, message_length):
        mdi3.write(str(I2[i])+"\n")
        mdq3.write(str(Q2[i])+"\n")
        
    I3 = ss.filtfilt(b, a, I2)
    Q3 = ss.filtfilt(b, a, Q2)
    
#   output I3 and Q3 to files for visualisation during testing
    for i in range(0, len(I3)):
        mdi4.write(str(I3[i])+"\n")
        mdq4.write(str(Q3[i])+"\n")
        
#   sync sample clock to bit clock
    clock_sync_list = I2[0:int((samples_per_symbol*4)-1)]
    largest_amplitude = biggest_magnitude(clock_sync_list)
    start_index = clock_sync_list.index(largest_amplitude)
    
    for centre_sample in range(start_index, message_length, int(samples_per_symbol)):
        #centre_sample is the sample at the centre of the integration range
        #slice_start is the first sample in integration range
        #slice_end is the end sample in the integration range
        slice_start = int(centre_sample - (samples_per_symbol/2))
        slice_end = int(centre_sample + (samples_per_symbol/2))
        iarm_samples = I2[slice_start:slice_end]
        qarm_samples = Q2[slice_start:slice_end]
        
        Isamples_sum = sum(iarm_samples)
        Qsamples_sum = sum(qarm_samples)
        
        if message_sync != True:
            
            if len(sync_list) == 8:
                del sync_list[0:2]
                
            #decide on bit for Q arm
            if Qsamples_sum > 0:
                sync_list.append("0")
            elif Qsamples_sum < 0:
                sync_list.append("1")
            else:
                sync_list.append("?")
            
            #decide on bit for I arm
            if Isamples_sum > 0:
                sync_list.append("1")
            elif Isamples_sum < 0:
                sync_list.append("0")
            else:
                sync_list.append("?")

        
            sync_word = ''.join(sync_list)[::]
            if sync_word == sync_key:
                message_sync = True
                print("SYNCED")
            else:
                print("not synced")
                
        #once output is synced with the start of a byte, start output to output file
        else:
            if Qsamples_sum > 0:
                output_file.write("0")
            elif Qsamples_sum < 0:
                output_file.write("1")
            else:
                output_file.write("?")
         
                 # I bit output
            if Isamples_sum > 0:
                output_file.write("1")
            elif Isamples_sum < 0:
                output_file.write("0")
            else:
                output_file.write("?")
 
    
#   close all open files
    encoded_message.close()
    output_file.close()
    mdi1.close()
    mdq1.close()
    mdi2.close()
    mdq2.close()
    mdi3.close()
    mdq3.close()
    mdibits.close()
    mdqbits.close()
    
    return 0


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

# =============================================================================
# For testing
# =============================================================================
def generate_timer_file(start=0, sample_rate=SAMPLE_RATE, encoded_message_file=ENCODED_MESSAGE):
    timer_file = open("timer.txt", "w")
    sample_interval = 1/sample_rate
    length_of_file = count_lines(encoded_message_file)
    
    t = start
    for i in range(0, length_of_file):
        timer_file.write(str(t)+"\n")
        t+=sample_interval
    
    timer_file.close()

    return 0

def sanity_checks():
    if CARRIER_FREQUENCY > SAMPLE_RATE:
        print("carrier_frequency greater than sample rate")
        return 1
    elif CARRIER_FREQUENCY < BIT_RATE:
        print("carrier_frequency less than bit rate")
        return 1
    elif BIT_RATE > SAMPLE_RATE:
        print("bit rate greater than sample rate")
        return 1
    
    return 0

def count_lines(file_name):
    file = open(file_name, "r")
    line_count = 0
    for line in file:
        line_count+=1
    file.close()
    return line_count
    
# =============================================================================
# GO    
# =============================================================================
checks = sanity_checks()
if checks != 1:
    main()
    generate_timer_file(-1/BIT_RATE)