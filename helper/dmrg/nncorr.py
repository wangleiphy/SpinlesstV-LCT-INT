import pyalps
import argparse
from numpy import reshape , zeros , array , where 
import sys 

def nncorr(x, y, L):
    res = zeros(L/2+1)
    for ij, ninj in zip(x, y):
        dist = min( abs(ij[0] -ij[1]), L - abs(ij[0] -ij[1])) # torus distance 
        shellsize = 1. if (dist == L/2) else 2.
        #print ij[0], ij[1], dist, shellsize 

        res[dist] += (ninj - 0.25)/shellsize 

    res = res/L
    res[0] = 0.25 

    return res 

parser = argparse.ArgumentParser()
parser.add_argument("-fileheader", default='chainL32APBCXL32_W1_N16_V1.0_SWEEPS16_M1600', help="fileheader")
parser.add_argument("-dir", default='../../data/dmrg/', help="dir")
args = parser.parse_args()

resfiles = pyalps.getResultFiles(dirname= args.dir, prefix=args.fileheader)
resfiles.sort()
print resfiles 

data = pyalps.loadEigenstateMeasurements(resfiles,'nncorr')
data = pyalps.flatten(data)
#print data 

#print data,len(data)

for d in data:
    V = d.props['V0']
    L = d.props['L']

    res = nncorr(d.x, d.y[0], L)
    for x, y in enumerate(res):
        print x, y 

