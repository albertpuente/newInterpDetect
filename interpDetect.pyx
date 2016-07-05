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
        InterpDetection(int rows, int cols) except +
        void detect(unsigned short* vm, long t0, long t1)

def interpDetect(filePath):

    # Read data from a .brw (HDF5) file
    rf = openHDF5file(filePath)
    nFrames, samplingRate, nRecCh, chIndices = getHDF5params(rf)
    
    # Duration of the recording in seconds
    nSec = nFrames / samplingRate

    # Start detection
    cdef InterpDetection * SpkD = new InterpDetection(64, 64)
    tInc = min(500000, nFrames)
    tCut = 20

    # vm is indexed as follows:
    #     vm[channel + tInc*nChannels] or vm[i*chRows + j + tInc*nChannels]
    cdef np.ndarray[unsigned short, mode = "c"] vm = np.zeros(nRecCh*tInc, dtype=ctypes.c_ushort)
    
    # Iterate
    for t0 in range(0L, nFrames - tCut, tInc - tCut):
        t1 = min(t0 + tInc, nFrames)

        vm = readHDF5(rf, t0, t1) 
        SpkD.detect(&vm[0], t0, t1)
