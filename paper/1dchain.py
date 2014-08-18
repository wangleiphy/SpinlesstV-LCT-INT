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
parser.add_argument("-f1", nargs='+', default="params", help="fileheaders")
parser.add_argument("-f2", nargs='+', default="params", help="fileheaders")

parser.add_argument("-logscale", action='store_true',  help="logscale")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

###################################################################
resultFiles = []
for fileheader in args.f1:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)
resultFiles = list(set(resultFiles))
resultFiles.sort()

print resultFiles 

res = pyalps.loadMeasurements(resultFiles, 'nncorr')
res = pyalps.flatten(res)
#pyalps.propsort(res,'V') 

print res 

fig = plt.figure(figsize = (8, 8))
ax1 = plt.subplot(211)

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

pyalps.plot.plot(res)

if args.logscale:
    plt.gca().set_yscale('log')

###################################################################
resultFiles = []
for fileheader in args.f2:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)
resultFiles = list(set(resultFiles))
resultFiles.sort()

res = pyalps.loadMeasurements(resultFiles, 'nncorr')
res = pyalps.flatten(res)
data = pyalps.loadMeasurements(resultFiles, ['S2'])
#print data 
res = pyalps.collectXY(data, x='NA', y='S2', foreach = ['V'])

ax2 = plt.subplot(212)
icolor = 0
lines = []
for d in res:
    V = d.props['V']
    d.props['line'] = 's'
    d.props['xlabel'] = '$N_A$'
    d.props['ylabel'] = '$S_2$'
    d.props['label'] = r'$V/t=%g$'%(V)
    d.props['color'] = colors[icolor]

    NA, S2 = loadtxt('../data/dmrg/chainL32APBCX_V'+str(V)+'_S2.dat', unpack = True, comments= '#', usecols= (0,1))
    line = plt.plot(NA+1, S2, '-', color = colors[icolor])
    lines.append(line)

    icolor = (icolor +1)%len(colors)

pyalps.plot.plot(res)
plt.xlim([1,31])

legends = ['$V/t=1$','$V/t=2$','$V/t=3$','$V/t=4$']
fig.legend(lines, legends, loc = (0.15, 0.6), shadow=True,fancybox=True)


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


