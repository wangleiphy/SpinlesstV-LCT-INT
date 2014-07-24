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


V, M2, IntE, Energy = loadtxt('../data/ed/V_M2_honeycombL3W3APBC.dat', unpack = True, comments= '#', usecols= (0,1,2,3))


icolor = 0
#for d in res1:
#    d.props['label'] = r'$\langle k \rangle$'
#    d.props['ylabel'] = ''
#    d.props['line'] = 'o'
#    d.props['color'] = colors[icolor]
#    icolor = (icolor +1)%len(colors)


for d in res2:

    d.props['line'] = 'o'
    d.props['ylabel'] = ''
    d.props['label'] = r'$\langle\hat{V}\rangle$'
    d.props['color'] = colors[icolor]
    plt.plot(V, IntE, '-',  c = colors[icolor] )
    icolor = (icolor +1)%len(colors)


for d in res3:
    #L = d.props['L']

    d.props['line'] = 'o'
    d.props['ylabel'] = ''
    d.props['label'] = r'$\langle\hat{H}\rangle$'
    d.props['color'] = colors[icolor]
    plt.plot(V, Energy, '-', c = colors[icolor] )
    icolor = (icolor +1)%len(colors)

for d in res4:
    #L = d.props['L']

    d.props['line'] = 'o'
    d.props['ylabel'] = ''
    d.props['xlabel'] = '$V/t$'
    d.props['label'] = '$\chi$'
    d.props['color'] = colors[icolor]
    plt.plot(V, M2, '-', c = colors[icolor])
    icolor = (icolor +1)%len(colors)


#pyalps.plot.plot(res1)
pyalps.plot.plot(res2)
pyalps.plot.plot(res3)
pyalps.plot.plot(res4)

plt.legend(loc='best')

#plt.xlim([1.28,1.42])
#plt.ylim([0,0.4])


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
