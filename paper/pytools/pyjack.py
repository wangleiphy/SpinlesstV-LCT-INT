# MC Jackknife analysis
# by Troyer Group, Institute for Theoretical Physics, ETH Zurich (C) 2012 - 2013

import numpy as np

def prepare_jacknife(data,nbins=None):
    bs = 1
    if nbins == None:
       nbins = len(data)/bs
    bins = [np.mean(data[:bs*nbins],axis=0)]
    for b in range(nbins):
        bins.append((bins[0]*nbins-np.mean(data[bs*b:bs*(b+1)],axis=0)) / (nbins-1))
    return np.array(bins)

# return jack mean, error
def evaluate_jackbins(jacks):
    jacks = np.asanyarray(jacks)
    m = len(jacks)-1
    err = np.sqrt( (m-1) * np.mean((jacks[1:]-jacks[0])**2,axis=0) )
    bias = (m-1) * (np.mean(jacks[1:],axis=0) - jacks[0])
    mean_ = jacks[0]- bias
    return mean_, err


