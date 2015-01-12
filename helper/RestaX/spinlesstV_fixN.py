from numpy import zeros  , log , pi , exp 
from bitops import btest, bget, bset, bclear, fermionParity, popcount

def buildH(Kmat, Ntot, V):
    '''
    -t*hopping +V*(n-0.5) (n-0.5)
    '''
    
    Nsites = Kmat.shape[0]
    Nstates = 1 << Nsites
    
    index2state = {}
    state2index = {}
    index = 0
    for state in range(Nstates):
        if popcount(state) == Ntot:
            index2state[index] = state 
            state2index[state] = index 
            index += 1 
    
    print index2state 
    print state2index 
    
    H = zeros((len(index2state), len(index2state)), float)

    # loop over all basis states
    for index, state in index2state.iteritems():
        #print bin(state)

        for s in range(Nsites):
            for j in range(Kmat.indptr[s], Kmat.indptr[s+1]):
                sj = Kmat.indices[j]  # a neighboring site sj 

                #interaction , if to avoid double counting 
                if sj > s: H[index, index] += V* (bget(state, s) -0.5) *  (bget(state, sj) -0.5)

                #hopping from site s to site sj 
                if btest(state, s):
                    parity1 = fermionParity(state, s)
                    state2 = bclear(state,s)
                    if not btest(state2,sj):
                        parity2 = parity1 ^ fermionParity(state2,sj)
                        state3=bset(state2,sj)

                        H[index, state2index[state3]] += Kmat.data[j] * (1-2*parity2)
    
    return H , index2state 


def calc_RestaX(Kmat, Ntot, V):
    '''
    compute renyi EE S2 
    '''
    
    Nsites = Kmat.shape[0]
    Hmat, index2state = buildH(Kmat, Ntot, V)

    w, v = eigh(Hmat)

    print w[0:4]/Nsites 
    Nstates = len(w)

    res = 0.0
    # loop over all basis states
    for index, state in index2state.iteritems():
        phase = 0.0
        for si in range(Nsites):
            phase += si * (bget(state, si) -0.5) 
    
        res +=  exp(1J*4.*pi/Nsites*phase) * v[index, 0]**2

    return res 


if __name__=='__main__':
    from numpy import array , zeros , linspace , dot , exp , diag , arange 
    import scipy.sparse as sps 
    from numpy.linalg import eigh 
    import sys 

    Thop = 1.0
    L = 8
    N = L/2 

    #Kinetic energy matrix 
    Kmat = zeros((L, L),float)

    #PBC 
    for s in range(L):   
        Kmat[s, (s+1)%L] = -Thop 
        Kmat[(s+1)%L, s] = -Thop 
    #OBC 
    #for s in range(L-1):   
    #    Kmat[s, s+1] = -Thop 
    #    Kmat[s+1, s] = -Thop 
    
    print Kmat 

    Kmat = sps.csr_matrix(Kmat)

    #for V in arange(0.0, 4.0, 0.2):
    if True:
        V = 4.0 
        X = calc_RestaX(Kmat, N, V)

        print V, X.real, X.imag 
