import os.path 
import subprocess
import re , sys 
from math import cos , sin , pi 

def writeParameterFile(fname,parms):
    """ This function writes a text input file for simple ALPS applications like DMFT
    
        The arguments are:
        
          filename: the name of the parameter file to be written
          parms: the parameter dict
    """
    f = file(fname,'w')
    for key in parms:
      value = parms[key]
      if type(value) == str:
        f.write(str(key)+' = "' + value + '"\n')
      else:
        f.write(str(key)+' = ' + str(value) + '\n')
    f.close()
    return fname

'''
write parameters for main 
'''

def params(lattice, L, W, V= 1.0, T = 0.1, Maxorder = 2048, Nscratch=10, SWEEPS=1000000, THERMALIZATION=100000, NBLOCKS = 15, STEPS_PER_BLOCK=10, folder='../data/', textoutput=0):
    
    key = lattice.replace(' ','') 
    key += 'L' + str(L)\
           +'W' + str(W)\
           +'V'+str(V)\
           +'T' + str(T)\
           +'MAXORDER'+ str(Maxorder)\
           +'Therm'+str(THERMALIZATION)\
           +'Sweeps'+str(SWEEPS) \
           +'Nscratch'+str(Nscratch)\
           +'NBLOCKS'+ str(NBLOCKS)\
           +'STEPSPERBLOCK'+str(STEPS_PER_BLOCK)


    inputname = '../jobs/'+ key +'.in'
    outputname = folder + key +'.dat'
    
    parms ={'LATTICE_LIBRARY' : "../input/mylattices.xml"
            # above we should not change 

            ,'LATTICE'  : lattice
            ,'filename' : outputname
            ,'textoutput' :textoutput 
            ,'L'  : L 
            ,'W'  : W

            ,'MAX_ORDER' : Maxorder
            ,'V' : V
            ,'TEMPERATURE' : T
            ,'RECALC_PERIOD' : Nscratch 
            ,'THERMALIZATION' : THERMALIZATION
            ,'SWEEPS' : SWEEPS 

            ,'NBLOCKS'  : NBLOCKS 
            ,'STEPS_PER_BLOCK'  : STEPS_PER_BLOCK
            }

    writeParameterFile(inputname, parms)

    return inputname 
