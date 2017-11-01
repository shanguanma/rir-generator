# This must be run as follows on WSL. 
# export LD_LIBRARY_PATH=$PWD
# python3 example_1.py

import pyrirgen

c = 340                          # Sound velocity (m/s)
fs = 16000                       # Sample frequency (samples/s)
r = [2, 1.5, 2]                  # Receiver position [x y z] (m)
s = [2, 3.5, 2]                  # Source position [x y z] (m)
L = [5, 4, 6]                    # Room dimensions [x y z] (m)
rt = 0.4                         # Reverberation time (s)
n = 2048                         # Number of samples

h = pyrirgen.generateRir(L, s, r, soundVelocity=c, fs=fs, reverbTime=rt, nSamples=n)

# Uncomment if you want to save the generated RIR.
#import numpy as np
#outfile = 'example_1'
#h = np.array(h)
#np.savez(outfile, filter=h)

