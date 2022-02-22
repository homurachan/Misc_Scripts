# Misc_Scripts
Some useful scripts I wrote.

calc_cc_bfactor_sjors_param.py calculates the envelop of chromatic aberration with given voltage, cc, delta_E. The rest parameters are from Nakane, Takanori, et al. "Single-particle cryo-EM at atomic resolution." Nature 587.7832 (2020): 152-156.

convert_motioncor2_log_to_aln.py just converts bunches of log files from motioncor2 to the .aln files. You can use these .aln file in motioncor2 with '-InAln your_folder'.

read_relion_run_data_write_BILD_file.py converts the refined relion starfile to a BILD file. It runs quite slow if your angle_step is small.

relion31_star_to_30.py converts relion3.1 starfiles (especiallly refined starfiles) to relion3.0 ver.
