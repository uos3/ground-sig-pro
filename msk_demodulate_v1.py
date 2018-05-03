# =============================================================================
# msk_demodulate_v1
#
# combines frequency finder and demodulator
# =============================================================================

import math as m
import numpy as np
import scipy.signal as ss
import soundfile as sf

BIT_RATE = 6000
OUTPUT_FILE = "demodulated_bits.bin"

# =============================================================================
# define all functions used to demodulate signal
# =============================================================================

def run_demodulator(input_file, bit_rate=BIT_RATE):
    """
    Reads in the amplitudes from the .wav file to a list, and then calls the 
    find_carrier() to retrieve the carrier function and msk_decode() to demodulate
    the message.
    
    Input
            input_file       -   .wav file of MSK modulated signal
            bit_rate         -   bit rate of incoming message
    
    Output
            does not return anythin
            creates .BIN file containint demodulated bits
    
    Variables
            signal_in       -   array containing MSK signal sample amplitudes
            sample_rate     -   sample rate
            fc              -   carrier frequency
            bit_rate        -   bit rate of incoming signal
            OUTPUT_FILE     -   filename to output demodulated bits to
    
    Functions called
            find_carrier()
            msk_decode()
            
            
    """
    signal_in, sample_rate = sf.read("msk_encoded_message.wav")
    
    fc = find_carrier(signal_in, sample_rate)
    betais = msk_decode(fc, sample_rate, bit_rate, signal_in, OUTPUT_FILE)
    
    return betais

def find_carrier(signal_in, sample_rate):
    """
    Takes in list of MSK signal sample amplitudes annd outputs carrier frequency
    
    Input
            signal_in       -   array containing MSK signal sample amplitudes
            sample_rate     -   sample rate
    
    Output
            fc              -   carrier frequency
    
    Variables
            n               -   number of signal samples
            sp              -   array of complex signal strengths 
            freq            -   array containing sample frequencies
            sp_abs          -   absolute value of sp
    
    Functions called
            len()
            numpy.fft.fft()
            numpy.fft.fftfreq()
            abs()
            numpy.argmax()            
    """
    # read in data from input array
    n = len(signal_in)
    
    # perform fast fourier transform
    sp = np.fft.fft(signal_in)
    # generate list of frequencies
    freq = np.fft.fftfreq(n, 1/sample_rate)
    # abosulte values of above lists
    sp_abs = abs(sp)
    freq = abs(freq)
    #find frequency correlating to maximum signal strength
    fc = freq[np.argmax(sp_abs)]

    return fc

def msk_decode(fc, sample_rate, bit_rate, data_in, output_file):
    """
    Takes carrier frequency and MSK signal sample amplitudes and performs
    non coherent demodulation to retrieve the modulated message.
    
    Input
            fc                  -   carrier frequency
            sample_rate         -   sample rate
            bit_rate            -   bit rate
            data_in             -   list containing MSK sample amplitudes
            output_file         -   name of .BIN file to output bits to
    
    Output
            does not return anythin
            creates .BIN file containint demodulated bits   
    
    Variables
            encoded_list        -   list containing input MSK signal ampitudes
            decoded_message     -   file object to store decoded message
            m_phases            -   number ofbits in one symbol
            sample_period       -   sample_period
            samples_per_symbol  -   samples per symbol (2 * samples per bit)
            bit_period          -   bit period
            samples_per_bit     -   samples per bit
            wc                  -   cyclic carrier frequency (2*pi*fc)
            wm                  -   cyclic modulation frequency (pi/(2*Tb))
            B_ORDER1            -   order of butterworth filter
            B_CUTOFF1           -   cutoff frequency of butterworth filter (ratio of cyclic modulation frequency and sample frequency)
            mdb                 -   numerator polynomial of butterworth filter
            mda                 -   denominator polynomial of butter filter
            
            R                   } - lists to hold signal samples after first   
            I                   }   multiplication in demodulator algorithm
            
            Rcs                 }
            Rss                 }
            Ics                 }
            Iss                 } - lists to hold signal samples after second
            RcTs                }   multiplication in demodulator algorithm
            RsTs                }
            IcTs                }
            IsTs                }
    
            ts                  -   list to hold demodulation clock timesteps
            t                   -   demodulation clock at a specific time
            message_length      -   number of samples in signal
            
            Rcss                }
            Rsss                }
            Icss                }
            Isss                } - lists to hold signal samples after low pass
            RcTss               }   filter 
            RsTss               }
            IcTss               }
            IsTss               }
            
            centre_sample       -   index in signal list that is the centre
                                    sample of current symbol being demodulated
            bit_1_start         -   index of signal list where current symbol
                                    being demodulated starts                                    
            bit_2_end           -   index of signal list where current symbol
                                    being demodulated ends
            
            Rc                  -   integration of Rcss over one bit period
            Rs                  -   integration of Rsss over one bit period
            Ic                  -   integration of Icss over one bit period
            Is                  -   integration of Isss over one bit period
            RcT                 -   integration of RcTss over one bit period  
            RsT                 -   integration of RcTss over one bit period 
            IcT                 -   integration of IcTss over one bit period
            IsT                 -   integration of IsTss over one bit period
                
            C1                  } - results of sum junctions
            D1                  } 
            C2                  -   C1 squared
            D2                  -   D1 squared
            
            beta1               } 
            beta2               } - symbol likelihoods
            beta3               }
            beta4               }
    
            output_symbol       -   variable to hold symbol to ouput to output file
            
            betais              -   list containing bit likelihoods
        
            beta_max            -   maximum bit likelihood
            max_index           -   index of betais corresponding to bet_max
        
            output_symbols      -   list of the four possible symbols to output
            output_symbol       -   the correct symbol to output for the current 
    
    Functions called
            list()
            open()
            int()
            len()
            range()
            sum()
            max()
            list.append()
            math.sin()
            math.cos()
            file.write()
            file.close()
            scipy.signal.butter()
            scipy.signal.filtfilt
            what_betais()
    """
#   read in data from input array and convert to list
    encoded_list = list(data_in)

#   open output file
    decoded_message = open(output_file, "w")
    
#   constants used throughout decoding process
    m_phases = 2
    sample_period = 1/sample_rate
    samples_per_symbol = int(m_phases*sample_rate/bit_rate)
    bit_period = 1/bit_rate
    samples_per_bit = int(sample_rate/bit_rate)
    wc = 2*m.pi*fc
    wm = m.pi/(2*bit_period)
    
    #butterworth filter constants
    B_ORDER1 = 5
    B_CUTOFF1 = (m.pi/(2/bit_rate))/sample_rate
   
#   define butterworth filter coefficients
    mdb, mda = ss.butter(B_ORDER1, B_CUTOFF1, 'low')
    
#   define variables to hold data during decoding process
    R = []
    I = []
    Rcs = []
    Rss = []
    Ics = []
    Iss = []
    RcTs = []
    RsTs = []
    IcTs = []
    IsTs = []
    ts = []
    beta1 = []
    beta2 = []
    beta3 = []
    beta4 = []
    output_symbol = ""
    
    t = 0
    ts = [0]
    message_length = len(encoded_list)
    
    for i in range(0, message_length):
        t += sample_period
        ts.append(t)
        
        #Perform first multiplications for R and I
        R.append(encoded_list[i]*m.cos(wc*ts[i]))
        I.append(encoded_list[i]*m.sin(wc*ts[i]))
        
        #Perform second multiplications for 8 variables
        Rcs.append(R[i]*m.cos(wm*(ts[i]-bit_period)))
        Rss.append(R[i]*m.sin(wm*(ts[i]-bit_period)))
        Ics.append(I[i]*m.cos(wm*(ts[i]-bit_period)))
        Iss.append(I[i]*m.sin(wm*(ts[i]-bit_period)))
        RcTs.append(R[i]*m.cos(wm*ts[i]))
        RsTs.append(R[i]*m.sin(wm*ts[i]))
        IcTs.append(I[i]*m.cos(wm*ts[i]))
        IsTs.append(I[i]*m.sin(wm*ts[i]))
    
    #Pass each channel through a low pass filter
    Rcss = ss.filtfilt(mdb, mda, Rcs)
    Rsss = ss.filtfilt(mdb, mda, Rss)
    Icss = ss.filtfilt(mdb, mda, Ics)
    Isss = ss.filtfilt(mdb, mda, Iss)
    RcTss = ss.filtfilt(mdb, mda, RcTs)
    RsTss = ss.filtfilt(mdb, mda, RsTs)
    IcTss = ss.filtfilt(mdb, mda, IcTs)
    IsTss = ss.filtfilt(mdb, mda, IsTs)
    
    #Demodulate
    for centre_sample in range(samples_per_bit, message_length, samples_per_symbol):
        bit_1_start = int(centre_sample - samples_per_bit)
        bit_2_end = int(centre_sample + samples_per_bit)      
   
        Rc = sum(Rcss[centre_sample:bit_2_end])
        Rs = sum(Rsss[centre_sample:bit_2_end])
        Ic = sum(Icss[centre_sample:bit_2_end])
        Is = sum(Isss[centre_sample:bit_2_end])
        RcT = sum(RcTss[bit_1_start:centre_sample])
        RsT = sum(RsTss[bit_1_start:centre_sample])
        IcT = sum(IcTss[bit_1_start:centre_sample])
        IsT = sum(IsTss[bit_1_start:centre_sample])
        
    # =============================================================================
    # Calculation of beta1
    # =============================================================================
        C1 = RcT + IsT - Rs + Ic
        D1 = RsT - IcT + Rc + Is
        C2 = C1**2
        D2 = D1**2
        beta1 = C2 + D2
     
    # =============================================================================
    # Calculation of beta2
    # =============================================================================
        C1 = RcT + IsT + Rs + Ic
        D1 = RsT - IcT + Rc - Is
        C2 = C1**2
        D2 = D1**2
        beta2 = C2 + D2
        
    # =============================================================================
    # Calculation of beta3
    # =============================================================================
        C1 = RcT - IsT + Rs - Ic
        D1 = -RsT - IcT - Rc - Is
        C2 = C1**2
        D2 = D1**2
        beta3 = C2 + D2
    
    # =============================================================================
    # Calculation of beta4 
    # =============================================================================
        C1 = RcT - IsT - Rs - Ic
        D1 = -RsT - IcT - Rc + Is
        C2 = C1**2
        D2 = D1**2
        beta4 = C2 + D2

        # Order beta values
        betais = what_betais(fc, output_symbol, beta1, beta2, beta3, beta4)
        
        beta_max = max(betais)
        max_index = betais.index(beta_max)
        
        output_symbols = ["01", "00", "10", "11"]
        output_symbol = output_symbols[max_index]
        decoded_message.write(output_symbol)

    #close all open files
    decoded_message.close()
    
    return betais

def what_betais(fc, last_symbol, beta1, beta2, beta3, beta4):
    """
    This function is a fix to correct the consisent wrong output of the beta array order.
    If fc is between 6623 and 7968 each bit needs to be flipped
    Likewise if the last bit decoded is a certain bit
    
    Input
        fc                  -   carrier frequency
        last_symbol         -   the last symbol output to the output file
        
        beta1               } 
        beta2               } - symbol likelihoods
        beta3               }
        beta4               }
    
    Output
        beta1               } 
        beta2               } - symbol likelihoods
        beta3               }
        beta4               }
    
    Variables
        f_flip              -   multiplier to change beta values order depending on carrier frequency
        b_flip              -   multiplier to change beta values order depending on last symbol
    """
    if fc >= 6623 and fc <=7698:
        f_flip = -1
        if last_symbol == "11" or last_symbol == "01":
            b_flip  = -1
        else:
            b_flip = 1
    else:
        f_flip = 1
        if last_symbol == "00" or last_symbol == "10":
            b_flip  = -1
        else:
            b_flip = 1
    
    if f_flip*b_flip == 1:
        return [beta1, beta2, beta3, beta4]
    else:
        return [beta3, beta4, beta1, beta2]
    
    return 0

# =============================================================================
# Run
# =============================================================================
print(run_demodulator("msk_encoded_message.wav"))