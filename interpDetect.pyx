# distutils: language = c++
# distutils: sources = SpkDslowFilter.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
from ctypes import CDLL
import ctypes
from readUtils import openHDF5file, getHDF5params, readHDF5

cdef extern from "SpkDslowFilter.h" namespace "SpkDslowFilter":
    cdef cppclass InterpDetection:
        InterpDetection(int rows, int cols, double samplingRate) except +
        void detect(unsigned short* vm, int t0, int t1)

def interpDetect(filePath):

    # Read data from a .brw (HDF5) file
    rf = openHDF5file(filePath)
    nFrames, samplingRate, nRecCh, chIndices = getHDF5params(rf)

    nFrames = 2000 # <------------------------------------------------------------- DELETE

    # Duration of the recording in seconds
    nSec = nFrames / samplingRate
    
    # Start detection
    cdef InterpDetection * SpkD = new InterpDetection(64, 64, samplingRate)
    tInc = min(100000, nFrames)
    tCut = 20

    # vm is indexed as follows:
    #     vm[channel + tInc*nChannels] or vm[i*chRows + j + tInc*nChannels]
    cdef np.ndarray[unsigned short, mode = "c"] vm = np.zeros(nRecCh*tInc, dtype=ctypes.c_ushort)
    
    # Iterate
    for t0 in range(0, nFrames - tCut, tInc - tCut):
        t1 = min(t0 + tInc, nFrames)
        print 'Detecting frames:', t0, 'to', t1, '...'
        vm = readHDF5(rf, t0, t1) 
        SpkD.detect(&vm[0], t0, t1)
    
