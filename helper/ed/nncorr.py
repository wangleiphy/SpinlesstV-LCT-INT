import pyalps
import argparse
from numpy import reshape , zeros , array , where 
from pltgraph import parse  
import scipy.sparse as sps 
import sys 


def M2(nnmat):
    nsite = nnmat.shape[0] 

    res =0.0
    for si in range(nsite):
        for sj in range(nsite):
            parity = (-1)**(si + sj)
            res += parity * nnmat[si, sj]

    return res/(nsite*nsite)

def IntE(Kmat, nnmat):
    nsite = nnmat.shape[0] 
    res = 0.0
    for si in range(nsite):
        for j in range(Kmat.indptr[si], Kmat.indptr[si+1]): # neighboring of si 
            sj = Kmat.indices[j]
            res += nnmat[si, sj]
    
    return res/nsite 
 

parser = argparse.ArgumentParser()
parser.add_argument("-fileheader", default='honeycombL3W3L3_W3_N9_V', help="fileheader")
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
    
    ###forming the connection graph, we need it to compute the interaction energy 
    graph = parse('honeycombL3W3.graph')

    #read lattice 
    Row = []
    Col = []
    Typ = []
    edges = graph[1]
    for edge in edges.values():
        Row.append(int(edge['source'])-1)
        Col.append(int(edge['target'])-1)
        Typ.append(int(edge['type']))

    Nsite = max(max(Row), max(Col)) + 1

    #Kinetic part 
    Val = zeros(len(Typ), float)
    Typ = array(Typ)
    Val[where(Typ==0)]= -1.0 # actually hopping 
    Val[where(Typ==1)]= -1.0 # actually hopping 

    Kmat = sps.csr_matrix((Val, (Row, Col)), shape=(Nsite, Nsite)) # we have one term per bond, Kmat is not hermitian now 
    ###forming the connection graph 

    nnmat = reshape(d.y, (2*L*W, 2*L*W)) # <n_i n_j>

    print V, M2(nnmat), V*IntE(Kmat, nnmat) -1.5*V*0.25, e.y[0]/(2.*L*W) - 1.5*V * 0.25 
