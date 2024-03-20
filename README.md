# RIR-Generator

The image method, proposed by Allen and Berkley in 1979 [1], is probably one of the most frequently used methods in the acoustic signal processing community to create synthetic room impulse responses. 

This repository includes both Cython and NumPy implementations of the image method along with a tutorial written by Dr. Emanuel A. P. Habets. (The NumPy version is 20-30x slower than the Cython verion. This was written for my personal learning purposes.)

This repository was forked from https://github.com/Marvin182/rir-generator. For more information, visit the original repository. 

1. J.B. Allen and D.A. Berkley, "Image method for efficiently simulating small-room acoustics," Journal Acoustic Society of America, 65(4), April 1979, p 943.


Installing pyrirgen (Cython version) on Windows
-------------
- Python 3.5: pip install .\pyrirgen-py35-win_amd64\pyrirgen-0.0.0-cp35-cp35m-win_amd64.whl
- Python 3.6: pip install .\pyrirgen-py36-win_amd64\pyrirgen-0.0.0-cp36-cp36m-win_amd64.whl




Installing nprirgen (NumPy version)
-------------
pip install -e ./nprirgen


Installing pyrirgen (py version)
pip install -e ./pyrirgen

This package is used by `https://github.com/jsalt2020-asrdiar/jsalt2020_simulate/blob/master/libaueffect/room_simulators/genrir.py#L191`

[UPDATE](2024-3-20):
https://github.com/Marvin182/rir-generator is not exist now.
