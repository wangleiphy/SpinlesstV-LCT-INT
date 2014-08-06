import subprocess 
import time 
import re 
from numpy import arange, array 

Add = 0.4 
Remove = 0.4 

Maxorder = 8192
BETA = 40. #projection time 
BCmodifier = '' #only APBCX will have effect anything else will not affet the lattice 

#latticename = 'chain lattice'
latticename = 'honeycomb lattice'
###############################
nickname = 'PBC1'

Llist = [3]
Wlist = Llist 
#Vlist = [1.0, 1.4]
Vlist = arange(1.3, 1.41, 0.01)
#Vlist = arange(0.2, 1.6, 0.2)

itime_max = 1<<31
RECALC_PERIOD = 17
WRAP_REFRESH_PERIOD = 25

STEPS_PER_BLOCK = 1
NBLOCKS = 1024
THERMALIZATION = 10**5
SWEEPS = 5*10**6 
MEASUREMENT_PERIOD = 13        # in unit of block

##############################
wtime = '01:00:00'
tmin = 60
tmax = 600
ncores = 320  # a multiply of ntasks_per_node 
prog = '../bin/main_halfTheta'

resfolder = '/mnt/lnec/lewang/spinlessctbssdata/' + nickname  + '/'
h, m, s = [int(i) for i in wtime.split(':')]
Tlimit = max(3600*h + 60*m + s - int(tmax*2) , 0)
prog += ' -i '+ str(tmin) + ' -a ' + str(tmax) + ' -T ' + str(Tlimit) # + ' -c '

def submitJob(bin,args,jobname,wtime,run=False,ncores=20, wait=None):

#SBATCH --ntasks=%g
    #prepare the job file 
    job='''#!/bin/bash -l
#
#SBATCH --exclusive
#SBATCH --nodes=%g
#SBATCH --time=%s
#SBATCH --partition=dphys_largemem 
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

