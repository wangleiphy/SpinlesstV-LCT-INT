from glob import glob 
from os import rename

for fname in glob('/scratch/rosa/lewang/spinlessctbssdata/zeroT_checkpoint/honeycomblatticeAPBCXL9W9*'):
    rename(fname, fname.replace('Sweeps1000000', 'Sweeps2000000'))


