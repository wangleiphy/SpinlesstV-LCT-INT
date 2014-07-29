from glob import glob 
from os import rename, mkdir 
import os.path 
from shutil import move

import argparse

#'/scratch/rosa/lewang/spinlessctbssdata/zeroT_checkpoint/honeycomblatticeAPBCXL9W9*
#checkpoint_files = '/mnt/lnec/lewang/spinlessctbssdata/zeroTfine_checkpoint/honeycomblatticeAPBCXL12W12V1.0ITIMEMAX2147483648BETA40.0NBLOCKS1024STEPSPERBLOCK4WRAP10RECALC10MAXORDER8192Therm20000Sweeps1000000Nskip13*'


parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", default="/mnt/lnec/lewang/spinlessctbss/data/test", help="fileheaders")
args = parser.parse_args()

for fname in glob(args.fileheaders+'*.h5'):

    #newfname = fname.replace('Sweeps1000000', 'Sweeps2000000')
    #newfname = newfname.replace('STEPSPERBLOCK4', 'STEPSPERBLOCK2')

    #rename(fname, newfname)
    #print fname,  newfname
    
    if 'clone' in fname:

        fpath = fname[0: fname.find('.clone')] + '.chkp'

        if (not os.path.isdir(fpath)):
            mkdir(fpath)
    
        move(fname, fpath+'/'+fname[fname.find('clone'):])

        print fname 
        print fpath+'/'+fname[fname.find('clone'):]






