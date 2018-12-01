################ 2FSK Modulator ###############
# Uses 2-frequency shift keying, or BFSK[binary frequency shift keying]
# according to https://www.dsprelated.com/showarticle/1016.php

import math as m
import numpy as n
import soundfile as sf
import matplotlib.pyplot as plt

message = "error_test_file_input.txt"
carrier_f=10000
bit_rate = 6000
bit_period = 1 / bit_rate
sample_rate=48000
sample_period = 1/sample_rate
deviation = bit_rate/4 #maximum deviation
output_file_name = "msk_encoded_message.wav"

message_file = open(message, "r")
message = [i for i in message_file.read()]

phaseshift = [(((carrier_f-deviation)/bit_rate)%1)*2*m.pi,
(((carrier_f+deviation)/bit_rate)%1)*2*m.pi]
totalphaseshift=[0]


for i in range(0, len(message)):
    totalphaseshift.append(
        (totalphaseshift[i]+phaseshift[int(message[i])]) % (2*m.pi))
    message[i]=(2*int(message[i]))-1

def fsk_wave(t, carrier_f, bit_period, message, deviation):
    bit_phase = (t % bit_period) #how far into the bit we are
    bit_number =  int((t - bit_phase) / bit_period)
    current_phase = totalphaseshift[bit_number]
    current_freq = carrier_f+deviation*message[bit_number]
    current_sig = m.cos(2*m.pi*bit_phase*current_freq+current_phase)
    return current_sig

total_time = len(message) * bit_period
sample_times = n.arange(0, total_time-sample_period, sample_period)
samples = []
phaseshifts=[]
for i in range(0, len(sample_times)):
    samples.append(fsk_wave(i*sample_period, carrier_f, bit_period, message, deviation))
#plt.plot(sample_times, samples, 'ro')
#plt.show()
sf.write(output_file_name, samples, sample_rate)
