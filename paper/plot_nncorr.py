'''
python load_corr.py -f /cluster/work/scr6/lewang/PQMCDATA/nncorrhofstadter3latticep1periodicL6W6N36U -s
'''
import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import sys ,os  
import subprocess 
import socket
import subprocess

import argparse
from config import * 
from numpy import loadtxt , abs 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")
parser.add_argument("-y", default="nncorr", help="observable")

parser.add_argument("-copydata", action='store_true',  help="copy data")
parser.add_argument("-logscale", action='store_true',  help="logscale")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)
resultFiles = list(set(resultFiles))
resultFiles.sort()

#filter resultFilies
#for f in list(resultFiles):
#    if ('L3' in f) or ('Nneighbors3' not in f) :
#    if ('Therm1000000S' not in f) :
#    if ('Sweeps10000000Skip1000' not in f):
#    if ('Skip20' not in f):
#    if ('Sweeps1000000' in f):
#     if ('V1.42' in f):
#        resultFiles.remove(f)

print resultFiles 

res = pyalps.loadMeasurements(resultFiles, args.y)
res = pyalps.flatten(res)
#pyalps.propsort(res,'V') 

print res 

icolor = 0
for d in res:
    V = d.props['V']
    d.props['xlabel'] = r'$r$'

    if args.logscale:
        d.props['ylabel'] = r'$(-1)^r C(r)$'
        d.y = abs(d.y)
    else:
        d.props['ylabel'] = r'$C(r)$'

    d.props['label'] =  r'$V/t=%g$'%(V)

    d.props['line'] = 'o'
    d.props['color'] = colors[icolor]
    
    #R, nncorr = loadtxt('../data/dmrg/nncorr_chainL32_L32_W1_N16_V'+str(V)+'.dat', unpack=True, usecols = (0,1))
    R, nncorr = loadtxt('../data/dmrg/nncorr_chainL32APBCXL32_W1_N16_V'+str(V)+'.dat', unpack=True, usecols = (0,1))
    if args.logscale:
        plt.plot(R, abs(nncorr), '-', c=colors[icolor])
    else:
        plt.plot(R, nncorr, '-', c=colors[icolor])

    icolor = (icolor+1)%len(colors)


if args.copydata:
    for resultFile in resultFiles:
        cmd = ['cp', resultFile, '../data/']
        subprocess.check_call(cmd)


print pyalps.plot.convertToText(res)

pyalps.plot.plot(res)
plt.legend(loc='lower left')

if args.logscale:
    plt.gca().set_yscale('log')


if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    recipient = "lewang@phys.ethz.ch"
    message = 'Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
    message += '\n' + pyalps.plot.convertToText(res)
    subject = 'Figure: ' + args.outname

    machinename = socket.gethostname()
    if 'brutus' in machinename or 'monch' in machinename:
        pyalps.sendmail(recipient    # email address of recipients 
                       , subject = subject 
                       , message = message 
                       , attachment= args.outname 
                       )
    else:
        cmd = ['sendmail.py', '-t', recipient+',', '-s', 'Automatic email message from ALPS. '+ subject , '-m', message, '-a', args.outname]
        subprocess.check_call(cmd)

    os.system('rm '+args.outname)


