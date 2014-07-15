import subprocess 
from numpy import arange 

Add  = 0.3
Remove = 0.3

ZtoW = 0.2
WtoZ = 0.2

latticename = 'open chain lattice'
#latticename = 'honeycomb lattice'
#latticename = 'open honeycomb lattice'
#latticename = 'cylindrical honeycomb lattice'
#latticename = 'square lattice'
#latticename = 'piflux lattice'
###############################
nickname = 'AttractiveHubbard'

Llist = [6]
Wlist = [1]

NA0list = [0]
NA1list = [3]

#NAstep = 4 

#Tlist = [0.1]
#Tlist = arange(0.6, 1.2, 0.2)
Tlist = [0.2, 0.4]
#Vlist = arange(0.1, 1.6, 0.2)
#Vlist = arange(0.5, 10., 0.5)
#Vlist = arange(1.7, 2.4, 0.1)
Ulist = -arange(1., 11., 1.)
Mulist = [1.0]

Ntau = 1000
NSKIP = 100 
THERMALIZATION = 10**5
SWEEPS = 10**6
Nscratch = 1000
##############################

tmin = 60
tmax = 300
ncores = 16 
wtime = '1:00'
bin = '../bin/HubbardRenyiQMC'

resfolder = '/cluster/work/scr6/lewang/renyidata/' + nickname  + '/'
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
