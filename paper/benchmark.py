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
from matplotlib.ticker import MaxNLocator
import re 
from config import * 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)

print resultFiles 

data = pyalps.loadMeasurements(resultFiles, ['PertOrder','M2', 'IntE', 'Energy'])
#print data 
#res = pyalps.ResultsToXY(data, 'PertOrder', 'Walltime', foreach = ['V']) 

res1 = pyalps.collectXY(data, x='V', y='PertOrder', foreach = ['L'])
res2 = pyalps.collectXY(data, x='V', y='IntE', foreach = ['L'])

res3 = pyalps.collectXY(data, x='V', y='Energy', foreach = ['L'])
res4 = pyalps.collectXY(data, x='V', y='M2', foreach = ['L'])

#print pyalps.plot.convertToText(res)


V, M2, IntE, Energy = loadtxt('../data/ed/V_M2_honeycombL3W3.dat', unpack = True, comments= '#', usecols= (0,1,2,3))


icolor = 0
#for d in res1:
#    d.props['label'] = r'$\langle k \rangle$'
#    d.props['ylabel'] = ''
#    d.props['line'] = 'o'
#    d.props['color'] = colors[icolor]
#    icolor = (icolor +1)%len(colors)


fig = plt.figure()
ax1 = fig.add_subplot(111)

for d in res2:
    d.props['line'] = 'o'
    d.props['ylabel'] = ''
    d.props['label'] = r'$\langle\hat{H}_1\rangle/N$'
    d.props['color'] = colors[icolor]
    plt.plot(V, IntE, '-',  c = colors[icolor] )
    icolor = (icolor +1)%len(colors)


for d in res3:
    #L = d.props['L']
    d.props['line'] = 's'
    d.props['xlabel'] = '$V/t$'
    d.props['ylabel'] = 'Energy per site'
    d.props['label'] = r'$\langle\hat{H}\rangle/N$'
    d.props['color'] = colors[icolor]
    plt.plot(V, Energy, '-', c = colors[icolor] )
    icolor = (icolor +1)%len(colors)

#pyalps.plot.plot(res1)
pyalps.plot.plot(res2)
pyalps.plot.plot(res3)
plt.legend(loc='upper left')

plt.ylim([-1.2,0.3])

ax2 = ax1.twinx()

for d in res4:
    #L = d.props['L']
    d.props['line'] = '^'
    d.props['ylabel'] = 'CDW Structure Factor'
    d.props['label'] = '$M_2$'
    d.props['color'] = colors[icolor]
    plt.plot(V, M2, '-', c = colors[icolor])
    icolor = (icolor +1)%len(colors)

pyalps.plot.plot(res4)

#for tl in ax2.get_yticklabels():
#    tl.set_color(colors[icolor-1])
#ax2.yaxis.label.set_color(colors[icolor-1])

plt.legend(loc='upper right')

plt.xlim([0.0, 2.0])

plt.subplots_adjust(right=0.88)

################inset###################
#inset = plt.axes([0.22, 0.42, 0.26, 0.22])
#data = pyalps.loadMeasurements(resultFiles, ['PertOrder','Walltime'])
#res = pyalps.ResultsToXY(data, 'PertOrder', 'Walltime') 

#for d in res:
#    #d.y /= d.y[0].mean
#    d.y /=3600. 
#    d.props['xlabel'] = r'$\langle k \rangle$'
#    d.props['line'] = '-o'
#    d.props['color'] = colors[icolor]
#    icolor = (icolor+1)%len(colors)

#pyalps.plot.plot(res)
#plt.xlabel( r'$\langle k \rangle$', fontsize=14)
#plt.ylabel(r'Wall time (h)', fontsize=14)
#inset.xaxis.set_major_locator(MaxNLocator(4))
#inset.yaxis.set_major_locator(MaxNLocator(4))
#plt.ylim([0, 8])
#plt.xlim([0, 360])
################inset###################


if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    recipient = "lewang@phys.ethz.ch"
    message = 'Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
    #message += '\n' + pyalps.plot.convertToText(res)
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
