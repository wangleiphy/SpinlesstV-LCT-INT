#!/usr/bin/env python
'''
python load_vertices.py -f ../data/test.clone0.h5
'''

import h5py 
from numpy import  array 
import argparse

import pyalps
import matplotlib.pyplot as plt 
import pyalps.plot
import os , sys

parser = argparse.ArgumentParser()
parser.add_argument("-file", help="h5file")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="All.pdf",  help="output pdf file")

args = parser.parse_args()

h5 = h5py.File(args.file,'r')

Theta = float(h5['/parameters/BETA'][()])
itime_max = float(h5['/parameters/ITIME_MAX'][()])

vt = array(h5['/simulation/realizations/0/clones/0/checkpoint/vt'][()])
ilist = array(h5['/simulation/realizations/0/clones/0/checkpoint/vi'][()])
jlist = array(h5['/simulation/realizations/0/clones/0/checkpoint/vj'][()])

tlist = vt/itime_max *Theta 
print tlist 

plt.scatter(tlist,ilist, color='b',marker=".")
plt.scatter(tlist,jlist, color='b',marker=".")

#for t, i, j in zip(tlist, ilist, jlist):
#    plt.plot([t, t], [i, j], c='b', lw=2)
    

h5.close()

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
