This folder contains config files used for the evaluations of the DNA-Aeon manuscript (currently in review).
The files can also serve as example configuration files for specific error rates.

The first Number of the filename (0.5, 1.5, 2.0 etc.) describes the multiplier that was applied to the 
error rates described using a high mutagenesis kit by Press et al. (10.1073/pnas.2004821117)

The base rates are as follows:
Subsitutions: 0.0238
Deletions:    0.0082
Insertions:   0.0039

I.e. configs with a 2.0 as prefix were evaluated with a total error rate of:

2.0\*0.238 + 2.0\*0.0082 + 2.0\*0.0039 = 0.0718

IMPORTANT: The evaluations were carried out on a powerful computer (500 cores, 2 TB Ram). If you want to use the 
config files, you will probably have to adjust the *threads* parameter (and maybe the *size*) parameter. The threads parameter
has to be adjusted for both your available CPU threads and for your available memory.
