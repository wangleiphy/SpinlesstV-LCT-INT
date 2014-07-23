import pyalps
import argparse

def Asublattice(nnlist):
    # a list only contains correlation between the A site in origin cell with other sites 
    l = []
    for i in range(len(nnlist)):
        if (i%4==0) or (i%4==1):
            l.append(nnlist[i])

    return l 

def M2(l):
   
    res =0.0
    for s, ninj in enumerate(l):
        res += (-1)**s * ninj 

    return res/len(l)

def IntE(l):# since we might not have rotation symmetry, we add three nearest neighbors 
    return l[1] + l[2*W-1] + l[4*W-1]
 

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
    W = int(d.props['W'])
    #print d.x 
    #print d.y 
    l =  Asublattice(d.y[0])
    print V, M2(l), 0.5*V*IntE(l) -1.5*V*0.25, e.y[0]/(2.*L*W) - 1.5*V * 0.25 
    #for x, y in zip(d.x, d.y):
    #    print  '(', x, ') : ', y
