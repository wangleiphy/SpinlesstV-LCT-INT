import pyalps
from subprocess import check_call 
from math import sqrt , cos , sin , pi 
from numpy import arange 

L = 3
W = 3
N = L*W
Thop = 1.0

parms = []
for V in arange(0.0, 2.2, 0.01):
    parms.append({
          'LATTICE_LIBRARY'           : 'mylattices.xml',
          'MODEL_LIBRARY'             : 'mymodels.xml',
# 	      'LATTICE'                   : "honeycomb lattice", 
          'LATTICE'                   : 'honeycombL3W3', 
          'MODEL'                     : "spinless fermions",
	      'CONSERVED_QUANTUMNUMBERS'  : 'N',
  	      'TRANSLATION_SYMMETRY'      : 'false',
          'L'                         : L,
          'W'                         : W,
          't0'                        : Thop,  
          't1'                        : Thop,  
          'V0'                        : V,
          'V1'                        : V,
	      'N_total'                   : N, 
          'NUMBER_EIGENVALUES'        : 1,
          'BCx'                       : "periodic",
          'BCy'                       : "periodic", 
#          'MEASURE_CORRELATIONS[corr]': 'cdag:c',
          'MEASURE_CORRELATIONS[nncorr]': 'n:n'
#          'MEASURE_LOCAL[Nloc]'       : 'n',
#          'INITIAL_SITE'              : 0
          #"PRINT_EIGENVECTORS"        : 1
        })


input_file = pyalps.writeInputFiles('spinlesstV',parms)
#res = pyalps.runApplication('sparsediag',input_file,writexml=False)

folder = '../../data/ed/'
for p in parms:
    parmname = folder + str(p['LATTICE']).replace(" ", "")+'L'+str(p['L'])+'_W'+str(p['W'])+'_N'+str(p['N_total'])+'_V'+str(p['V0']) 

    input_file = pyalps.writeInputFiles(parmname, [p])
    #pyalps.runApplication('sparsediag',input_file) #,writexml=True)#,MPI=2)
    #check_call(['echo', 'bsub', "-R", '"rusage[mem=1024]"', '-oo',input_file.replace('.in.xml','.log'),'-W','1:00','sparsediag',input_file])
    check_call(['bsub', '-oo',input_file.replace('.in.xml','.log'),'-W','1:00','sparsediag',input_file])
    #check_call(['bsub', '-n', '5','-oo',input_file.replace('.in.xml','.log'),'-W','08:00','mpirun', 'sparsediag', '--mpi', input_file])
