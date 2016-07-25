#import pyximport
#pyximport.install()

from interpDetect import interpDetect

rawpath = '/disk/scratch/apuente/data/'
fileName = '30.brw'

interpDetect(rawpath + fileName)