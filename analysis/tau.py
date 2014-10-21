'''
python\ tau.py\ -f\ /cluster/work/scr6/lewang/PQMCDATA/Qpipiskip/squarelatticeperiodicL4W4N8U1Theta20.0dtau0.05Ntau1000Nscratch10Therm20000Sweeps1000000Skip
'''

import pyalps
#import pyalps.dwa 
#import matplotlib.pyplot as plt
#import pyalps.plot
from report import report

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")
args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)

resultFiles = list(set(resultFiles))

for filename in resultFiles:

    obslist = ['PertOrder', 'M2', 'IntE', 'Energy','nncorr','dlogKinEdV','dlogM2dV']

    print filename 
    #print pyalps.dwa.tau(filename, obslist)
    #print pyalps.dwa.thermalized(filename, obslist)
    report(filename, obslist)
    print "#######################################################################"
