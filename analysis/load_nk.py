#!/usr/bin/env python
'''
python load_nk.py -f nk.h5 -o nk.pdf
'''

from plotnk import plot_nk 
import h5py 
from numpy import imag , array 
import argparse

import sys

parser = argparse.ArgumentParser()
parser.add_argument("-file", help="h5file")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="All.pdf",  help="output pdf file")

args = parser.parse_args()

h5 = h5py.File(args.file,'r')
#kxlist = array(h5['Params']['kxlist'])

nk = array(h5['nk'][()])
print nk.shape 
h5.close()

if args.show:
    plot_nk(nk)
else:
    plot_nk(nk, output=args.outname) 
