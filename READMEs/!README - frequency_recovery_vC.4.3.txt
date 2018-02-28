Frequency Recovery Loop Readme

Current version vC.4.3
- does not work
- currently working on outputting quantifiable error from the costas_loop() function 
- produces extra files that will not be in final product for visualisation during development

#####
2 Main Functions
#####


1. frequency_finder

Inputs 				- VCO_start 		(starting gues for carrier frequency)
				- encoded_message 	(msk modulated message file)
				- sample rate

Outputs 			- VCO_f 		(carrier frequency)

Generated files 		- none



2. costas_loop

Inputs 				- VCO_f 		(guess for carrier frequency)
				- filt1a
				- filt1b
				- filt2a
				- filt2b
				- sample_rate
				- max_samples 		(amount of samples to inlcude when finding freqency)
				- bit period

Outputs 			- VCO_f 		(guess for carrier frequency that was input)
				- f_err 		(error)

Generated files 		- frvC_A1.txt		)
    				- frvC_B1.txt		)
    				- frvC_C1.txt		)
    				- frvC_A2.txt		)
    				- frvC_B2.txt		)
    				- frvC_A3.txt		) - files to visualise the signal at different points during frequency recovery
    				- frvC_B3.txt		)
    				- frvC_C2.txt		)
    				- frvC_C3.txt		)
    				- frvC_C4.txt		)
    				- frvC_timer.txt	)
