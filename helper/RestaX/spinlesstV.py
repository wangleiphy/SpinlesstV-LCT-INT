from numpy import zeros  , log , pi , exp 
from bitops import btest, bget, bset, bclear, fermionParity

def buildH(Kmat, V):
    '''
    -t*hopping +V*(n-0.5) (n-0.5)
    '''
    
    Nsites = Kmat.shape[0]
    Nstates = 1 << Nsites
    H = zeros((Nstates, Nstates), float)

    # loop over all basis states
    for state in range(Nstates):
        #print bin(state)

        for s in range(Nsites):
            for j in range(Kmat.indptr[s], Kmat.indptr[s+1]):
                sj = Kmat.indices[j]  # a neighboring site sj 

                #interaction , if to avoid double counting 
                if sj > s: H[state, state] += V* (bget(state, s) -0.5) *  (bget(state, sj) -0.5)

                #hopping from site s to site sj 
                if btest(state, s):
                    parity1 = fermionParity(state, s)
                    state2 = bclear(state,s)
                    if not btest(state2,sj):
                        parity2 = parity1 ^ fermionParity(state2,sj)
                        state3=bset(state2,sj)

                        H[state, state3] += Kmat.data[j] * (1-2*parity2)
    
    return H 

def calc_RestaX(Kmat, V, beta):
    '''
    compute renyi EE S2 
    '''
    
    Nsites = Kmat.shape[0]
    Hmat = buildH(Kmat, V)

    w, v = eigh(Hmat)
        
    w -= w[0]
    weights = exp(-beta * w)
    Z = weights.sum()  
    weights = weights/Z 

    Nstates = len(w)

    res = 0.0
    # loop over all basis states
    for n in range(Nstates):
        for state in range(Nstates):
            phase = 0.0
            for si in range(Nsites):
                phase += si * bget(state, si)
        
            res +=  exp(1J*2.*pi/Nsites*phase) * v[state, n]**2 * weights[n]

    return res 


def calc_RestaX(Kmat, V):
    '''
    compute renyi EE S2 
    '''
    
    Nsites = Kmat.shape[0]
    Hmat = buildH(Kmat, V)

    w, v = eigh(Hmat)

    print w[0:4]/Nsites 
    Nstates = len(w)

    res = 0.0
    # loop over all basis states
    for state in range(Nstates):
        phase = 0.0
        for si in range(Nsites):
            phase += si * (bget(state, si) -0.5) 
    
        res +=  exp(1J*2.*pi/Nsites*phase) * v[state, 0]**2

    return res 



if __name__=='__main__':
    from numpy import array , zeros , linspace , dot , exp , diag , arange 
    import scipy.sparse as sps 
    from numpy.linalg import eigh 
    import sys 

    Thop = 1.0
    L = 8

    #Kinetic energy matrix 
    Kmat = zeros((L, L),float)

    #PBC 
    #for s in range(L):   
    #    Kmat[s, (s+1)%L] = -Thop 
    #    Kmat[(s+1)%L, s] = -Thop 
    #OBC 
    for s in range(L-1):   
        Kmat[s, s+1] = -Thop 
        Kmat[s+1, s] = -Thop 
    
    print Kmat 

    Kmat = sps.csr_matrix(Kmat)

    for V in arange(0.0, 4.0, 0.2):
        X = calc_RestaX(Kmat, V)

        print V, X.real, X.imag 
