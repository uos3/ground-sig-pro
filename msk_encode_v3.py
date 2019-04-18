################ 2FSK Modulator ###############
# Modulates a binary signal using continuous-phase binary
#   frequency shift keying (CP-BFSK) aka Minimum Shift Keying (MSK),
#   according to the formulas at
#   https://www.dsprelated.com/showarticle/1016.php

import math as m
import numpy as n
import soundfile as sf
import matplotlib.pyplot as plt

message_filename = "encoder_input.txt"
carrier_frequency = 10000

bit_rate = 6000
bit_period = 1 / bit_rate

sample_rate = 480000
sample_period = 1/sample_rate

#The minimum deviation required from carrier frequency, as given in the article.
deviation = bit_rate/4
output_filename = "encoder_output.wav"

message_file = open(message_filename, "r")
message = [i for i in message_file.read()]


#The "phaseshift" array stores the difference in signal phase induced by
#   both a '0' bit and a '1' bit.

phaseshift = [(((carrier_f-deviation)/bit_rate)%1)*2*m.pi,
(((carrier_f+deviation)/bit_rate)%1)*2*m.pi]


#The "totalphaseshift" array keeps track of the phase of the signal
# at the beginning of each bit. The loop also converts the message from
# the [0, 1] binary form to the [-1, 1] form.

totalphaseshift=[0]
for i in range(0, len(message)):
    totalphaseshift.append(
        (totalphaseshift[i]+phaseshift[int(message[i])]) % (2*m.pi))
    message[i]=(2*int(message[i]))-1


#The "fsk_wave" function calculates the relative signal (from 0 to 1) at
#   a sampled time.
def fsk_wave(t):
    bit_number = int(t // bit_period) #which bit we're on
    bit_phase = t % bit_period #how far into the bit we are

    current_freq = carrier_f+deviation*message[bit_number] #freq at this bit
    
    beginning_phase = totalphaseshift[bit_number] #phase at beginning of bit
    current_phase = 2*m.pi*bit_phase*current_freq+beginning_phase
        #phase at signal point
    current_signal = m.cos(current_phase)
    return current_signal

total_time = len(message) * bit_period
sample_times = n.arange(0, total_time-sample_period, sample_period)
samples = []
phaseshifts=[]
for i in range(0, len(sample_times)):
    samples.append(fsk_wave(i*sample_period))
plt.plot(sample_times, samples, 'ro')
plt.show()
sf.write(output_filename, samples, sample_rate)
