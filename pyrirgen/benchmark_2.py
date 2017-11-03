# This must be run as follows on WSL. 
# export LD_LIBRARY_PATH=/mnt/c/Users/tayoshio/Documents/GitHub/rir-generator/python/
# python3 example_1.py

import pyrirgen
import nprirgen
from timeit import default_timer as timer

import numpy as np

c = 340       # sound velocity (m/s)
fs = 16000    # sample frequency (samples/s)
trials = 30  # number of trials

roomdim_range_x = [5, 10]     # width
roomdim_range_y = [5, 10]     # depth
roomdim_range_z = [2.5, 4.5]  # height

micarray = np.concatenate([np.zeros((1,3)), np.array([0.0425 * np.array([np.cos(i * np.pi/3), np.sin(i * np.pi/3), 0]) for i in range(6)])])

roomcenter_mic_dist_max_x = 0.5
roomcenter_mic_dist_max_y = 0.5
micpos_range_z = [0.6, 0.9]

spkr_mic_dist_range_x = [0.5, 4]
spkr_mic_dist_range_y = [0.5, 4]
spkrpos_range_z_rel2mic = [0.1, 0.5]

t60_range = [0.1, 0.4]

for i in range(trials):
    # Room dimensions
    L = np.array([np.random.uniform(*roomdim_range_x), 
                  np.random.uniform(*roomdim_range_y), 
                  np.random.uniform(*roomdim_range_z)])

    # Microphone position
    room_center = L / 2
    r = np.array([np.random.uniform(room_center[0] - roomcenter_mic_dist_max_x, room_center[0] + roomcenter_mic_dist_max_x),
                  np.random.uniform(room_center[1] - roomcenter_mic_dist_max_y, room_center[1] + roomcenter_mic_dist_max_y),
                  np.random.uniform(*micpos_range_z)])
    r = np.maximum(r, 0)
    r = np.minimum(r, L)
    R = micarray + r

    # Source position
    angle = np.random.uniform(0, 2 * np.pi)
    s = np.array([r[0] + np.random.uniform(*spkr_mic_dist_range_x) * np.cos(angle), 
                  r[1] + np.random.uniform(*spkr_mic_dist_range_y) * np.sin(angle), 
                  r[2] + np.random.uniform(*spkrpos_range_z_rel2mic)])
    s = np.maximum(s, 0)
    s = np.minimum(s, L)

    dist = np.linalg.norm(r - s)

    rt = np.random.uniform(*t60_range)  # T60
    n = int(rt * fs)  # number of samples

    # validity check
    V = np.prod(L)
    S = 2 * (L[0]*L[2] + L[1]*L[2] + L[0]*L[1])
    alpha = 24 * V * np.log(10) / (c * S * rt)
    if alpha > 1:
        print('[Trial #{:03}] Reflection coeffs were not able to be calculated from T60.'.format(i))
        continue


    # pyrirgen
    start = timer()
    h0 = pyrirgen.generateRir(L, s, R, soundVelocity=c, fs=fs, reverbTime=rt, nSamples=n)
    h0 = np.array(h0)
    end = timer()
    t_pyrirgen = end - start
   
    # nprirgen
    start = timer()
    h1, _, refl_coeffs = nprirgen.np_generateRir(L, s, R, soundVelocity=c, fs=fs, reverbTime=rt, nSamples=n)    
    end = timer()
    t_nprirgen = end - start

    # accuracy
    snr = 10 * np.log10(np.sum(np.abs(h0)**2) / np.sum(np.abs(h0 - h1)**2))

    print('[Trial #{:03}] pyrirgen: {:6.3f} sec \t nprirgen (rel): {:6.3f}\t snr={:6.3f} dB, L=[{:4.2f}, {:4.2f}, {:4.2f}], r=[{:4.2f}, {:4.2f}, {:4.2f}], s=[{:4.2f}, {:4.2f}, {:4.2f}], dis={:5.3f}, rt={:5.3f}, n={}, refl_coeffs=[{:4.2f}, {:4.2f}, {:4.2f}, {:4.2f}, {:4.2f}, {:4.2f}]'.format(i, t_pyrirgen, t_nprirgen/t_pyrirgen, snr, L[0], L[1], L[2], r[0], r[1], r[2], s[0], s[1], s[2], dist, rt, n, refl_coeffs[0], refl_coeffs[1], refl_coeffs[2], refl_coeffs[3], refl_coeffs[4], refl_coeffs[5]))

    