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

en = pyalps.loadEigenstateMeasurements(resfiles,'Energy')
en = pyalps.flatten(en)

#print data,len(data)

print '#V, M2, interaction energy per site, energy per site '
for d, e in zip(data, en):
    V = d.props['V0']
    L = d.props['L']
    W = d.props['W']
    print V, M2(d.y[0][::2]), 1.5*V*(d.y[0][1] - 0.25), e.y[0]/(2.*L*W) - 1.5*V * 0.25 
    #for x, y in zip(d.x, d.y):
    #    print  '(', x, ') : ', y
