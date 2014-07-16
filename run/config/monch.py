import subprocess 
import time 
import re 
from numpy import arange, array 

#Divisors[2^32 - 1] = {1, 3, 5, 15, 17, 51, 85, 255, 257, 771, 1285, 3855, 4369, 13107, 21845, 65535, 65537, 196611, 327685, 983055, 1114129, 3342387, 5570645, 16711935, 16843009, 50529027, 84215045, 252645135, 286331153, 858993459, 1431655765, 4294967295}

latticename = 'honeycomb lattice'
###############################
nickname = 'stablization'

Llist = array([6]) 
Tlist = 0.75/Llist 

Vlist = [1.3, 1.4]

RECALC_PERIOD = 10
UPDATE_REFRESH_PERIOD = 10
WRAP_REFRESH_PERIOD = 10

STEPS_PER_BLOCK = 1
NBLOCKS = 256
THERMALIZATION = 10**4
SWEEPS = 10**6 

##############################
wtime = '12:00:00'
tmin = 60
tmax = 600
ncores = 400  # a multiply of ntasks_per_node 
prog = '../bin/main'

resfolder = '/mnt/lnec/lewang/spinlessctbssdata/' + nickname  + '/'
h, m, s = [int(i) for i in wtime.split(':')]
Tlimit = max(3600*h + 60*m + s - int(tmax*2) , 0)
prog += ' -i '+ str(tmin) + ' -a ' + str(tmax) + ' -T ' + str(Tlimit) 

def submitJob(bin,args,jobname,wtime,run=False,ncores=20, wait=None):

#SBATCH --ntasks=%g
    #prepare the job file 
    job='''#!/bin/bash -l
#
#SBATCH --exclusive
#SBATCH --nodes=%g
#SBATCH --time=%s
#SBATCH --partition=dphys_compute
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-socket=10
#SBATCH --cpus-per-task=1
#SBATCH --job-name=%s
#SBATCH --output=%s
#SBATCH --error=%s'''%(ncores/20, wtime,jobname,jobname+'.log',jobname+'.log')

    if wait is not None:
        dependency ='''
#SBATCH --dependency=afterany:%d\n'''%(wait)
        job += dependency 


    job += '''
echo "The current job ID is $SLURM_JOB_ID"
echo "Running on $SLURM_JOB_NUM_NODES nodes:"
echo $SLURM_JOB_NODELIST
echo "Using $SLURM_NTASKS_PER_NODE tasks per node"
echo "A total of $SLURM_NTASKS tasks is used"\n'''

    job +='mpirun -rmk slurm '+ str(bin) + ' '
    for key, val in args.items():
        job += str(key) +' '+ str(val) + ' '

    #print job
    jobfile = open("jobfile", "w")
    jobfile.write("%s"%job)
    jobfile.close()

    #submit the job 
    if run:
        cmd = ['sbatch', 'jobfile']

        ret = subprocess.check_output(cmd)

        jobid = int(re.search(r'\d+', ret).group())
        print jobid , 'submitted'
        time.sleep(0.1)
        return jobid 

    else:

        cmd = ['cat','jobfile']
        subprocess.check_call(cmd)
        return None

