#import pyximport
#pyximport.install()

from interpDetect import interpDetect

rawpath = '/disk/scratch/apuente/data/'
fileName = '10.brw'

interpDetect(rawpath + fileName)