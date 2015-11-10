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

def params(lattice, BCmodifier, L, W, V= 1.0, BETA = 40., WINDOWSIZE= 4.0, Maxorder = 2048, itime_max=1073741824, MEASUREMENT_PERIOD = 10, RECALC_PERIOD=10,  WRAP_REFRESH_PERIOD=10, SWEEPS=1000000, THERMALIZATION=100000, NBLOCKS = 15, STEPS_PER_BLOCK=10, Add = 0.4, Remove = 0.4, folder='../data/', textoutput=0):
    
    key = lattice.replace(' ','') + BCmodifier  
    key += 'L' + str(L)\
           +'W' + str(W)\
           +'V'+str(V)\
           +'ITIMEMAX'+ str(itime_max)\
           +'BETA' + str(BETA)\
           +'WINDOWSIZE' + str(WINDOWSIZE)\
           +'NBLOCKS'+ str(NBLOCKS)\
           +'STEPSPERBLOCK'+str(STEPS_PER_BLOCK)\
           +'WRAP'+str(WRAP_REFRESH_PERIOD)\
           +'RECALC'+str(RECALC_PERIOD)\
           +'MAXORDER'+ str(Maxorder)\
           +'Therm'+str(THERMALIZATION)\
           +'Sweeps'+str(SWEEPS) \
           +'Nskip' + str(MEASUREMENT_PERIOD)

    totprob = Add + Remove
    if (totprob >=1.0 ):
        print 'Add + Remove= ', totalprob , ' are you sure ?'
        sys.exit(1)


    inputname = '../jobs/'+ key +'.in'
    outputname = folder + key +'.dat'
    
    parms ={'LATTICE_LIBRARY' : "../input/mylattices.xml"
            # above we should not change 

            ,'LATTICE'  : lattice
            ,'filename' : outputname
            ,'textoutput' :textoutput 

            ,'BCmodifier' : BCmodifier
            ,'L'  : L 
            ,'W'  : W

            ,'ITIME_MAX' : itime_max
            ,'MAX_ORDER' : Maxorder
            ,'V' : V
            ,'BETA' : BETA
            ,'WINDOWSIZE' : WINDOWSIZE

            ,'THERMALIZATION' : THERMALIZATION
            ,'SWEEPS' : SWEEPS 
            ,'MEASUREMENT_PERIOD' : MEASUREMENT_PERIOD

            ,'RECALC_PERIOD' : RECALC_PERIOD
            ,'WRAP_REFRESH_PERIOD' : WRAP_REFRESH_PERIOD

            ,'NBLOCKS'  : NBLOCKS 
            ,'STEPS_PER_BLOCK'  : STEPS_PER_BLOCK
            ,'Add'    : Add
            ,'Remove' : Remove  
            ,'MEASURE_M4' : 1
            }

    writeParameterFile(inputname, parms)

    return inputname 
