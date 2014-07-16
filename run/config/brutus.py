import subprocess 
from numpy import arange, array 

latticename = 'honeycomb lattice'
###############################
nickname = 'firsttry'

Llist = array([6]) 
Tlist = 0.75/Llist 

Vlist = [1.3, 1.4]

RECALC_PERIOD = 20
UPDATE_REFRESH_PERIOD = 20
WRAP_REFRESH_PERIOD = 50

STEPS_PER_BLOCK = 5
NBLOCKS = 51 
THERMALIZATION = 10**3
SWEEPS = 10**5
##############################

tmin = 60
tmax = 300
ncores = 16 
wtime = '1:00'
bin = '../bin/main'

resfolder = '/cluster/work/scr6/lewang/spinlessctbssdata/' + nickname  + '/'
#h, m = [int(i) for i in wtime.split(':')]
#Tlimit = max(3600*h + 60*m - int(tmax*2) , 0)

prog = 'mpirun '+ bin  + ' -i '+ str(tmin) + ' -a ' + str(tmax) 

def submitJob(bin,args,jobname,wtime,run=False,ncores=None, wait=[]):

        if run:
            cmd = ['bsub']
        else:
            cmd = ['echo', 'bsub']

        if ncores != None:
            cmd += ['-n', str(ncores)]

        cmd += ['-W',str(wtime)] 
        cmd += ['-J', jobname]  # LSF jobname
        cmd += ['-oo', jobname+'.log'] #log file 
        
        if wait is not None:
            if len(wait) > 0:
                conds='ended('+wait[0]+')'
                for w in wait[1:]: conds += '&&ended('+w+')'
                cmd += ['-w',conds]

        #cmd += ['-R', '\"rusage[mem=4096]\"'] #memory usage 

        cmd += [bin] 

        for key, val in args.items():
            cmd += [key, val]

        subprocess.check_call(cmd)
        #time.sleep(2)
        return None 
