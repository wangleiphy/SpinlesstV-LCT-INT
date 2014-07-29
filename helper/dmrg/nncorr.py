import pyalps
import argparse
from numpy import reshape , zeros , array , where 
import sys 

def nncorr(x, y, L):
    res = zeros(L/2+1)
    res[0] = 0.25 
    for ij, ninj in zip(x, y):
        dist = abs(ij[0] -ij[1])
        shellsize = 2. if (dist == L/2) else 1.
        res[dist] += (ninj - 0.25)/shellsize 

    return res/L 

parser = argparse.ArgumentParser()
parser.add_argument("-fileheader", default='openchainlatticeL32_W1_N16_V1.0_SWEEPS16_M800', help="fileheader")
parser.add_argument("-dir", default='../../data/dmrg/', help="dir")
args = parser.parse_args()

resfiles = pyalps.getResultFiles(dirname= args.dir, prefix=args.fileheader)
resfiles.sort()
print resfiles 

data = pyalps.loadEigenstateMeasurements(resfiles,'nncorr')
data = pyalps.flatten(data)
print data 

#print data,len(data)

for d in data:
    V = d.props['V0']
    L = d.props['L']

    print nncorr(d.x, d.y, L)

