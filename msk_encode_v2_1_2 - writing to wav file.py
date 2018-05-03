# =============================================================================
# MSK Encode v2_1_2
# Takes inpur file and modulates onto carrier wave using MSK modulation.
# Outputs a .wav file of modulated signal
# ==============================

import math as m
import soundfile as sf

# =============================================================================
# Modulator function
# =============================================================================
def msk_encode(carrier_A, carrier_f, bit_rate, sample_rate, message, output_file_name):
    message_file = open(message, "r")

    m_phases = 2
    message = [i for i in message_file.read()]
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
    sf.write(output_file_name, sum_channel, 48000)

    #close all files
    message_file.close()

    return 0

# =============================================================================
# Generates timer file for plotting output signal
# =============================================================================
def generate_timer_file(start, sample_rate, encoded_message_file):
    timer_file = open("timer.txt", "w")
    sample_interval = 1/sample_rate
    length_of_file = count_lines(encoded_message_file)
    
    t = start
    for i in range(0, length_of_file):
        timer_file.write(str(t)+"\n")
        t+=sample_interval
    
    timer_file.close()

    return 0    

def count_lines(file_name):
    file = open(file_name, "r")
    line_count = 0
    for line in file:
        line_count+=1
    file.close()
    return line_count
# =============================================================================
# run modulator
# =============================================================================
msk_encode(1, 10000, 6000, 48000, "error_test_file_input.txt", "msk_encoded_message.wav")
#generate_timer_file(0, 48000, "msk_encoded_message.txt")
print("msk_encode_v1 imported")