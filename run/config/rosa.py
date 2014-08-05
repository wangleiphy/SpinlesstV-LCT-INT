import subprocess 
import time 
from numpy import arange 

Add = 0.4 
Remove = 0.4 

Maxorder = 8192
BETA = 40. #projection time 
BCmodifier = '' #only APBCX will have effect anything else will not affet the lattice 

latticename = 'honeycomb lattice'
###############################
nickname = 'PBC'

Llist = [15]
Wlist = Llist 
#Vlist = [1.4]
Vlist = arange(1.3, 1.4, 0.01)
#Vlist = arange(0.2, 1.6, 0.2)

itime_max = 1<<31
RECALC_PERIOD = 11
WRAP_REFRESH_PERIOD = 25

STEPS_PER_BLOCK = 2
NBLOCKS = 1024
THERMALIZATION = 10**5
SWEEPS = 5*10**6  
MEASUREMENT_PERIOD = 13        # in unit of block

wtime = '24:00:00'
tmin = 300
tmax = 600
ncores = 320
prog = '../bin/main'
#######################################

resfolder = '/scratch/rosa/lewang/spinlessctbssdata/' + nickname  + '/'
h, m, s = [int(i) for i in wtime.split(':')]

Tlimit = max(3600*h + 60*m + s - int(tmax*2.) , 0)
prog += ' -i '+ str(tmin) + ' -a ' + str(tmax) + ' -T ' + str(Tlimit) + ' -c '

def submitJob(bin,args,jobname,wtime,run=False,ncores=None, wait = None):

#SBATCH --ntasks-per-node=16
#SBATCH --mem=2048
            #prepare the job file 
            job='''#!/bin/bash
#SBATCH --ntasks=%g
#SBATCH --time=%s
#SBATCH --account=s395
#SBATCH --job-name=%s
#SBATCH --output=%s
#SBATCH --error=%s\n'''%(ncores,wtime,jobname,jobname+'.log',jobname+'.log')

            job +='aprun -n '+  str(ncores)+' '+ str(bin) + ' '
            for key, val in args.items():
                job += str(key) +' '+ str(val) + ' '
            
            #print job
            jobfile = open("jobfile", "w")
            jobfile.write("%s"%job)
            jobfile.close()
            
            #submit the job 
            if run:
                cmd = ['sbatch','jobfile']
            else:
                cmd = ['cat','jobfile']

            subprocess.check_call(cmd)
            time.sleep(0.05)

