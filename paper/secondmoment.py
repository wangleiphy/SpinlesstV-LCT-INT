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
from numpy import array , sqrt 
from config import * 
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from NRlist import NRlist 


import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")
parser.add_argument("-copydata", action='store_true',  help="copy data")

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
    L = int(d.props['L'])
    #d.props['xlabel'] = r'$R$'
    #d.props['label'] =  r'$V=%g$'%(V)

    r = pyalps.DataSet()
    up  = 0.
    down = 0.
    CRmax = abs(d.y[-1])
    for x, y in enumerate(d.y):
        up += (y* (-1)**x-CRmax)* x*x #* NRlist[L][x]
        down += (y* (-1)**x-CRmax)    #* NRlist[L][x]

    r.y = array([sqrt(up/down)/L])

    r.props = d.props 
    r.props['observable'] = 'xioverL'

    res.append(r)

print res 
res = pyalps.collectXY(res, x='V', y='xioverL', foreach = ['L'])
pyalps.propsort(res,'L')

#scale data 
icolor = 0
for d in res:
    L = d.props['L']
    d.y = array([y*L**args.b for y in d.y])

    d.props['xlabel'] = r'$V$'
    d.props['ylabel'] = r'$\xi/L$'
    d.props['label'] = '$L=%g$' %(L)
    d.props['line'] = '-o'
    d.props['color'] = colors[icolor]
    icolor = (icolor+1)%len(colors)

if args.copydata:
    for resultFile in resultFiles:
        cmd = ['cp', resultFile, '../data/']
        subprocess.check_call(cmd)

#pyalps.propsort(res,'V') 
print pyalps.plot.convertToText(res)

fig = plt.figure(figsize = (8, 8))

ax1 = plt.subplot(211)
pyalps.plot.plot(res)
ax1.yaxis.set_major_locator(MaxNLocator(4))
#plt.axvline(args.xc,color='k')
#ax1.get_xaxis().set_visible(False)
plt.legend(loc='best')

#########################
at = AnchoredText("a",prop=dict(size=18), frameon=True,loc=1,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

#further scale x 
for d in res:
    L = d.props['L']
    d.x = (d.x- args.xc) * L ** args.a
    d.props['xlabel'] = r'$(V-V_c)L^{1/\nu}$'

ax2 = plt.subplot(212)
pyalps.plot.plot(res)
ax2.yaxis.set_major_locator(MaxNLocator(4))
plt.legend(loc='best')
#ax2.yaxis.set_major_locator(MaxNLocator(4))

#########################
at = AnchoredText("b",prop=dict(size=18), frameon=True,loc=1,)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.gca().add_artist(at)
#########################

plt.subplots_adjust(hspace=0.2)


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


