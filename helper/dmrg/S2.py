import pyalps
import argparse
from numpy import reshape , zeros , array , where 
import sys 


parser = argparse.ArgumentParser()
parser.add_argument("-fileheader", default='chainL32APBCXL32_W1_N16_V1.0_SWEEPS16_M1600', help="fileheader")
parser.add_argument("-dir", default='../../data/dmrg/', help="dir")
args = parser.parse_args()

resfiles = pyalps.getResultFiles(dirname= args.dir, prefix=args.fileheader)
resfiles.sort()
print resfiles 

data = pyalps.loadEigenstateMeasurements(resfiles,'Renyi2')
data = pyalps.flatten(data)

#print data 

for d in data:
    #print d.y 
    for i , S2 in enumerate(d.y[0]):
        print i, S2 
