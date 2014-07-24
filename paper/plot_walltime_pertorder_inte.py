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

#filter resultFilies
#for f in list(resultFiles):
#    L = int(re.search('L([0-9]*)W',f).group(1)) 
#    W = int(re.search('W([0-9]*)N',f).group(1)) 
#    N = int(re.search('N([0-9]*)U',f).group(1)) 
#    V= float(re.search('V([0-9]*\.?[0-9]*)ITIMEMAX',f).group(1)) 
#    if (N!= L*W):  
#        resultFiles.remove(f)
    
#    if V< 1.3 or V>1.4:
#        resultFiles.remove(f)


print resultFiles 

data = pyalps.loadMeasurements(resultFiles, ['PertOrder','Walltime', 'IntE'])
#print data 
#res = pyalps.ResultsToXY(data, 'PertOrder', 'Walltime', foreach = ['V']) 

res1 = pyalps.collectXY(data, x='V', y='PertOrder', foreach = ['L'])
res2 = pyalps.collectXY(data, x='V', y='IntE', foreach = ['L'])

res3 = pyalps.collectXY(data, x='V', y='Walltime', foreach = ['L'])

#print pyalps.plot.convertToText(res)

icolor = 0
for d in res1:
    d.props['label'] = r'$\langle k \rangle$'
    d.props['ylabel'] = ''
    d.props['line'] = '-o'
    d.props['color'] = colors[icolor]
    icolor = (icolor +1)%len(colors)

for d in res2:
    L = d.props['L']
    W = d.props['W']
    Theta = d.props['BETA']

    d.y *= -(2*L*W)*Theta 
    d.props['line'] = '-o'
    d.props['label'] = r'$-\Theta\langle\hat{V}\rangle$'
    d.props['ylabel'] = ''
    d.props['color'] = colors[icolor]
    icolor = (icolor +1)%len(colors)

########################################################
for d in res3:
    #L = d.props['L']

    d.props['label'] = 'Wall time'
    d.props['ylabel'] = 'Wall time (s)'
    d.props['line'] = '-o'
    d.props['color'] = colors[icolor]
    icolor = (icolor +1)%len(colors)



fig = plt.figure()
ax1 = fig.add_subplot(111)

pyalps.plot.plot(res1)
pyalps.plot.plot(res2)

V, IntE = loadtxt('../data/ed/V_M2_honeycombL3W3APBC.dat', unpack = True, comments= '#', usecols= (0,2))
plt.plot(V, -IntE * 18. * 40.)

plt.legend(loc='upper left')

ax2 = ax1.twinx()
pyalps.plot.plot(res3)

plt.legend(loc='lower right')

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
