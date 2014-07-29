from glob import glob 
from os import rename

#'/scratch/rosa/lewang/spinlessctbssdata/zeroT_checkpoint/honeycomblatticeAPBCXL9W9*

checkpoint_files = '/mnt/lnec/lewang/spinlessctbssdata/zeroTfine_checkpoint/honeycomblatticeAPBCXL12W12V1.0ITIMEMAX2147483648BETA40.0NBLOCKS1024STEPSPERBLOCK4WRAP10RECALC10MAXORDER8192Therm20000Sweeps1000000Nskip13*'

for fname in glob(checkpoint_files):
    newfname = fname.replace('Sweeps1000000', 'Sweeps2000000')
    newfname = newfname.replace('STEPSPERBLOCK4', 'STEPSPERBLOCK2')

    rename(fname, newfname)
    #print fname,  newfname


