'''
python load_results.py -f ../data/pscanLhofstadter3latticep1periodicL*Ntau2000 -s 
'''
import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import os, sys 
import subprocess 
import socket
import argparse
from numpy import array , linspace , sqrt , arange , loadtxt 
from pytools.floatwitherror import FloatWithError as fwe
import re 
from config import * 
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from extrapolate import extrapolate

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

parser.add_argument("-extrapolate", action='store_true',  help="do curve fitting")
parser.add_argument("-nextrapolate", type=int, default= 4 ,  help="number of points used in extrapolation")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)
resultFiles = list(set(resultFiles))

#filter resultFilies
for f in list(resultFiles):
    L = int(re.search('L([0-9]*)W',f).group(1)) 
    V= float(re.search('V([0-9]*\.?[0-9]*)ITIMEMAX',f).group(1)) 
    
#    if V< 1.3 or V>1.4:
#        resultFiles.remove(f)

#    if L in [15]:
#        resultFiles.remove(f)
    #if L in [3]:
#        resultFiles.remove(f)


PBCfiles = []
ABPCfiles = []
for f in resultFiles:
    if 'APBCX' in f:
        ABPCfiles.append(f)
    else:
        PBCfiles.append(f)

data = pyalps.loadMeasurements(PBCfiles, 'Energy')
PBCres = pyalps.collectXY(data, x='L', y='Energy', foreach = ['V'])
pyalps.propsort(PBCres,'V')

data = pyalps.loadMeasurements(ABPCfiles, 'Energy')
APBCres = pyalps.collectXY(data, x='L', y='Energy', foreach = ['V'])
pyalps.propsort(APBCres,'V')

#print pyalps.plot.convertToText(res)

res1 = []
res2 = []
#scale data 
for d1, d2 in zip(PBCres, APBCres):

    V = d1.props['V']

    d1.x = 1./d1.x 
    d2.x = 1./d2.x 

    d1.props['label'] = 'QMC-PBC'
    d2.props['label'] = 'QMC-APBC'

    d2.props['xlabel'] = '$1/L\,\mathrm{or}\,1/D$'
    d1.props['ylabel'] = ''
    d2.props['ylabel'] = ''

    d1.props['color'] = colors[0]
    d2.props['color'] = colors[1]

    d1.props['line'] = '-o'
    d2.props['line'] = '-s'
 
    if V == 1.0:
        res1.append(d1)
        res1.append(d2)
    else:
        res2.append(d1)
        res2.append(d2)

fig = plt.figure(figsize = (8, 8))

ax1 = plt.subplot(211)

pyalps.plot.plot(res1)

D, En = loadtxt('../data/iPEPS/V1.0.dat', unpack=True, usecols = (0,1))
plt.plot(1./D, En, '-', marker = '*', c=colors[2], label = 'iPEPS', markersize=8)

if args.extrapolate:
    extrapolate(res1, args.nextrapolate)

ax1.yaxis.set_major_locator(MaxNLocator(5))
ax1.get_xaxis().set_visible(False)
plt.xlim([0, 0.25])
plt.ylim([-0.92, -0.88])

plt.title('$V/t=1.0$')
plt.legend(loc='upper left')

#########################
at = AnchoredText("a",prop=dict(size=18), frameon=True,loc=1,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

ax2 = plt.subplot(212, sharex=ax1)
pyalps.plot.plot(res2)
D, En = loadtxt('../data/iPEPS/V1.4.dat', unpack=True, usecols = (0,1))
plt.plot(1./D, En, '-', marker = '*', c=colors[2], label = 'iPEPS', markersize=8)

if args.extrapolate:
    extrapolate(res2, args.nextrapolate)

ax2.yaxis.set_major_locator(MaxNLocator(4))

plt.xlim([0, 0.25])
plt.ylim([-0.98, -0.95])

#########################
at = AnchoredText("b",prop=dict(size=18), frameon=True,loc=1,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

plt.title('$V/t=1.4$')
plt.legend(loc='upper left')

##shared y label 
yyl=plt.ylabel(r'Energy per site')
yyl.set_position((yyl.get_position()[0],1)) # This says use the top of the bottom axis as the reference point.
yyl.set_verticalalignment('center') 

plt.subplots_adjust(hspace =0.2,left=0.15)

if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    recipient = "lewang@phys.ethz.ch"
    message = 'Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
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
