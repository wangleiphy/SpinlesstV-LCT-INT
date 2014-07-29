import pyalps
from subprocess import check_call 
from math import sqrt , cos , sin , pi 
from numpy import arange 

L = 64 
W = 1
N = L*W/2 

Thop = 1.0
M = 3200

parms = []
#for V in arange(0.0, 2.2, 0.05):
for V in [1.0, 2.0, 3.0, 4.0]:
    parms.append({
          'LATTICE_LIBRARY'           : 'mylattices.xml',
          'MODEL_LIBRARY'             : 'mymodels.xml',
          'LATTICE'                   : 'chainL64APBCX', 
          'MODEL'                     : "spinless fermions",
	      'CONSERVED_QUANTUMNUMBERS'  : 'N',
  	      'TRANSLATION_SYMMETRY'      : 'false',
          'L'                         : L,
          'W'                         : W,
          't0'                        : Thop,  
          't1'                        : -Thop,  
          'V0'                        : V,
          'V1'                        : V,
	      'N_total'                   : N, 
          'NUMBER_EIGENVALUES'        : 1,
#          'MEASURE_CORRELATIONS[corr]': 'cdag:c',
          'MEASURE_CORRELATIONS[nncorr]': 'n:n', 
#          'MEASURE_LOCAL[Nloc]'       : 'n',
#          'INITIAL_SITE'              : 0
          #"PRINT_EIGENVECTORS"        : 1
          'SWEEPS'                     : 20, 
          'MAXSTATES'                  : M
        })


input_file = pyalps.writeInputFiles('spinlesstV',parms)

folder = '../../data/dmrg/'
for p in parms:
    parmname = folder + str(p['LATTICE']).replace(" ", "")+'L'+str(p['L'])+'_W'+str(p['W'])+'_N'+str(p['N_total'])+'_V'+str(p['V0']) + '_SWEEPS' + str(p['SWEEPS']) + '_M'  + str(p['MAXSTATES'])

    input_file = pyalps.writeInputFiles(parmname, [p])
    #pyalps.runApplication('mps_optim',input_file) #,writexml=True)#,MPI=2)

    #check_call(['echo', 'bsub', "-R", '"rusage[mem=1024]"', '-oo',input_file.replace('.in.xml','.log'),'-W','1:00','sparsediag',input_file])

    check_call(['bsub', '-oo',input_file.replace('.in.xml','.log'),'-W','24:00','mps_optim',input_file])

    #check_call(['bsub', '-n', '5','-oo',input_file.replace('.in.xml','.log'),'-W','08:00','mpirun', 'sparsediag', '--mpi', input_file])
