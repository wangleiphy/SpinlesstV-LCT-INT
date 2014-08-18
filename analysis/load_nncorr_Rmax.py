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
import re 
from numpy import array 

import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

parser.add_argument("-copydata", action='store_true',  help="copy data")
parser.add_argument("-logscale", action='store_true',  help="logscale")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

group = parser.add_argument_group()
group.add_argument("-xc", type = float,  default=0., help="xc")
group.add_argument("-a", type = float,  default=0., help="a")
group.add_argument("-b", type = float,  default=0., help="b")

args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)
resultFiles = list(set(resultFiles))

resultFiles.sort()

#filter resultFilies
for f in list(resultFiles):
    L = int(re.search('L([0-9]*)W',f).group(1)) 
    V= float(re.search('V([0-9]*\.?[0-9]*)ITIMEMAX',f).group(1)) 
    
#    if V not in [1.33, 1.34, 1.35, 1.36]:
#        resultFiles.remove(f)
    if V< 1.3 or V>1.4 or L in [3]:
        resultFiles.remove(f)


print resultFiles 

data = pyalps.loadMeasurements(resultFiles, 'nncorr')
data = pyalps.flatten(data)

res = []
for d in data:
    V = d.props['V']
    L = d.props['L']
    d.props['xlabel'] = r'$R$'
    d.props['label'] =  r'$V=%g$'%(V)

    
    r = pyalps.DataSet()
    r.y = [d.y[-1]]
    r.props = d.props 
    r.props['observable'] = 'farthestnncorr'

    res.append(r)

print res 
res = pyalps.collectXY(res, x='V', y='farthestnncorr', foreach = ['L'])

#for d in res:
#    d.y *= d.x**2
#    d.x = 1./d.x 

#scale data 
for d in res:
    L = d.props['L']
    d.x = (d.x- args.xc) * L ** args.a 
    d.y = array([y*L**args.b for y in d.y])

    d.props['xlabel'] = r'$(V-V_c)L^{a}$'
    d.props['ylabel'] = r'$C(R_\mathrm{max})L^{b}$'
    d.props['label'] = '$L=%g$' %(L)
    d.props['line'] = '-o'

if args.copydata:
    for resultFile in resultFiles:
        cmd = ['cp', resultFile, '../data/']
        subprocess.check_call(cmd)

#pyalps.propsort(res,'V') 
print pyalps.plot.convertToText(res)

pyalps.plot.plot(res)
plt.legend(loc='upper left')

#plt.xlim([0, 0.2])
#plt.ylim([0, 0.015])

if args.logscale:
    plt.gca().set_yscale('log')


#plt.figure()
#pyalps.plot.plot(farthest)
#plt.legend(loc='upper left')


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


