#!/usr/bin/env python
'''
python load_clone.py  -f /mnt/lnec/lewang/HofstadterHubbarddata/shift/squarelatticep1q4L4W4Npart8nshift0.48U-3.3ITIMEMAX2147483648BETA40.0NBLOCKS1024STEPSPERBLOCK1WRAP31RECALC41MAXORDER32768Therm500000Sweeps1000000Nskip17Add0.4Remove0.4.chkp/clone*0.h5 -x PertOrder -y Dbocc -s
'''

import h5py 
from numpy import  array 
import argparse

import pyalps
import matplotlib.pyplot as plt 
import pyalps.plot
import os , sys

parser = argparse.ArgumentParser()
parser.add_argument("-files", nargs = '+', help="h5file")

parser.add_argument("-x", default='PertOrder', help="x")
parser.add_argument("-y", default='k', help="y")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="All.pdf",  help="output pdf file")

args = parser.parse_args()

xlist = []
ylist = []
for f in args.files:
    h5 = h5py.File(f,'r')

    x = float(h5['/simulation/realizations/0/clones/0/measurements/'+args.x+'/mean/value'][()])
    y = float(h5['/simulation/realizations/0/clones/0/measurements/'+args.y+'/mean/value'][()])

    xlist.append(x)
    ylist.append(y)

    print x,y 

    h5.close()

plt.scatter(xlist, ylist)


if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    pyalps.sendmail('lewang@phys.ethz.ch'    # email address of recipients 
            , message='Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
                   , attachment= args.outname 
                   , subject='Figure: ' + args.outname 
                   )
    os.system('rm '+args.outname)
