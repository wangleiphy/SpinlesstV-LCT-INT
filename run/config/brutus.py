import subprocess 
from numpy import arange, array 

Add = 0.4
Remove = 0.4  

WINDOWSIZE = 4.0 # the window overwhich we perform measurment  
Maxorder = 16384 
NBLOCKlist = [1024, 2048, 4096, 8192]
BETA = 40. #projection time 
#BCmodifier = "APBCX"
BCmodifier = ""

#latticename = 'chain lattice'
latticename = 'square lattice'
#latticename = 'honeycomb lattice'
###############################
nickname = 'PBC'

Llist = [4, 8, 12, 16, 20]
Wlist = Llist 
#Vlist = [1.0, 2.0, 3.0, 4.0]
#Vlist = arange(1.2, 1.5, 0.02)
#Vlist = arange(0.2, 2.2, 0.2)
Vlist = arange(0.5, 5.0, 0.5)

itime_max = 1<<31
RECALC_PERIOD = 17
WRAP_REFRESH_PERIOD = 25 

STEPS_PER_BLOCK = 1
NBLOCKS = 1024
THERMALIZATION = 10**5            # in unit of block 
SWEEPS = 2*10**6                  # in unit of the the whole system 
MEASUREMENT_PERIOD = 33           # in unit of block
##############################

tmin = 60
tmax = 300
ncores = 16 
wtime = '18:00'
bin = '../bin/zeroT'

resfolder = '/cluster/work/scr6/lewang/spinlessctbssdata/' + nickname  + '/'
#h, m = [int(i) for i in wtime.split(':')]
#Tlimit = max(3600*h + 60*m - int(tmax*2) , 0)

prog = 'mpirun '+ bin  + ' -i '+ str(tmin) + ' -a ' + str(tmax) + ' -c '

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
        return [jobname]
