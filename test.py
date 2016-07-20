#import pyximport
#pyximport.install()

from interpDetect import interpDetect

rawpath = '/media/albert/OS/data/'
fileName = 'P29_16_05_14_retina02_left_stim3_fullarray_fullfieldHDF5.brw'

interpDetect(rawpath + fileName)