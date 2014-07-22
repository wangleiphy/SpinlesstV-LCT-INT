import pyalps
import argparse

def M2(ninjlist):
    res =0.0
    for s, ninj in enumerate(ninjlist):
        res += (-1)**s * ninj 

    return res/len(ninjlist)

   



parser = argparse.ArgumentParser()
parser.add_argument("-fileheader", default='honeycomblatticeL2_W2_N4_', help="fileheader")
parser.add_argument("-dir", default='../../data/ed/', help="dir")
args = parser.parse_args()

resfiles = pyalps.getResultFiles(dirname= args.dir, prefix=args.fileheader)
resfiles.sort()
print resfiles 

data = pyalps.loadEigenstateMeasurements(resfiles,'nncorr')
data = pyalps.flatten(data)

#print data,len(data)
for d in data:
    V = d.props['V0']
    print V, M2(d.y[0][::2]), 1.5*V*(d.y[0][1] - 0.25) # energy per site
    #for x, y in zip(d.x, d.y):
    #    print  '(', x, ') : ', y
