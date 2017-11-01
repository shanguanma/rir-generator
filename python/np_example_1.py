import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import nprirgen
import numpy as np
import matplotlib.pyplot as plt

c = 340                          # Sound velocity (m/s)
fs = 16000                       # Sample frequency (samples/s)
r = [2, 1.5, 2]                  # Receiver position [x y z] (m)
s = [2, 3.5, 2]                  # Source position [x y z] (m)
L = [5, 4, 6]                    # Room dimensions [x y z] (m)
rt = 0.4                         # Reverberation time (s)
n = 2048                         # Number of samples

h, _, _ = nprirgen.np_generateRir(L, s, r, soundVelocity=c, fs=fs, reverbTime=rt, nSamples=n)

plt.figure()
plt.plot(h)

# Uncomment to compare the result with that generated by pyrirgen.
#h0 = np.load('example_1.npz')
#h0 = h0['filter']
#plt.figure()
#plt.plot(h0 - h)

plt.show()
