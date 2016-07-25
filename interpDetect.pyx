# distutils: language = c++
# distutils: sources = SpkDslowFilter.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
from ctypes import CDLL
import ctypes
from readUtils import openHDF5file, getHDF5params, readHDF5

import os
clear = lambda: os.system('clear')
class bcolors:
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    

cdef extern from "SpkDslowFilter.h" namespace "SpkDslowFilter":
    cdef cppclass InterpDetection:
        InterpDetection(int rows, int cols, double samplingRate) except +
        void detect(unsigned short* vm, int t0, int t1, int tCut)

def interpDetect(filePath):

    # Read data from a .brw (HDF5) file
    rf = openHDF5file(filePath)
    nFrames, samplingRate, nRecCh, chIndices = getHDF5params(rf)

    # Duration of the recording in seconds
    nSec = nFrames / samplingRate
    
    # Start detection
    cdef InterpDetection * SpkD = new InterpDetection(64, 64, samplingRate)
    tInc = min(10000, nFrames)
    tCut = 50

    # vm is indexed as follows:
    #     vm[channel + tInc*nChannels] or vm[i*chRows + j + tInc*nChannels]
    cdef np.ndarray[unsigned short, mode = "c"] vm = np.zeros(nRecCh*tInc, dtype=ctypes.c_ushort)
    
    # Iterate
    
    for t0 in range(0, nFrames - tCut, tInc - tCut):
        t1 = min(t0 + tInc, nFrames)

        if os.name == 'posix': 
            displayProgress(t0, t1, nFrames)    

        
        print '\nDetecting frames:', t0, 'to', t1, ':'
        vm = readHDF5(rf, t0, t1) 
        SpkD.detect(&vm[0], t0, t1, tCut)


    print '\nDone.'
    
def displayProgress(t0, t1, nFrames):
    os.system('clear')
    col = int(os.popen('stty size', 'r').read().split()[1])
    g = np.divide(t0*col, nFrames)
    b = max(np.divide(t1*col, nFrames) - g, col - g)
    k = col - g - b
    print bcolors.GREEN + '█'*g + bcolors.BLUE +  '█'*b+ bcolors.ENDC +  '█'*k + '\n'