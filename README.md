README.md

## msk_encode_v2_1_2
Reads in .txt file nd outputs .wav file with MSK encoded message

carrier_A 		use 1

carrier_f 		use anything between 6000 and 14000. Expected carrier frequency will be 10000

bit_rate 		will be 6000

sample_rate		needs to be 48000

message 		.txt file name with message to encode: must be in the directory 

output_file_name	.wav file name for output MSK signal


### To run:
1. Make sure there is a message file in the same directory or specify the path in the script when you call the function.
2. Run the function, it will generate a .wav file.


## msk_demodulate_v1
Reads in the .wav file and ouputs demodulated bits

input_file 		this is the output from msk_encode_v2_1_2

bit_rate 		default set to 6000 as this is what the bit rate should be


It will generate a .BIN file with the demodulated message.

Notes: needs to operate at 48000 Hz sample rate, 6000 bps and the start of the message must be the sync key '00101010'

### To run:
1. Make sure the .wav file is in the same directory or specify the path in the script when you call the function
2. Run the function, it will generate a .BIN file
