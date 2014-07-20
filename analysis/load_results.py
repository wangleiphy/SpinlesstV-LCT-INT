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
from numpy import array , linspace , sqrt , arange 
import re 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

parser.add_argument("-x", default="L", help="variable")
parser.add_argument("-y", default="M2", help="observable")
parser.add_argument("-copydata", action='store_true',  help="copy data")


group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

group = parser.add_argument_group()
group.add_argument("-b", type = float,  default=0., help="b")


args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)

#filter resultFilies
#for f in list(resultFiles):
#    L = int(re.search('L([0-9]*)W',f).group(1)) 
#    W = int(re.search('W([0-9]*)N',f).group(1)) 
#    N = int(re.search('N([0-9]*)U',f).group(1)) 
#    V= float(re.search('latticeV([0-9]*\.?[0-9]*)T',f).group(1)) 
#    if (N!= L*W):  
#        resultFiles.remove(f)

#    elif ('L3' in f) or ('L9' in f) or ('L15' in f):
#    elif ('Theta155' in f):
#        resultFiles.remove(f)
#    if V not in [ 1.24,  1.26,  1.28,  1.3 ,  1.32,  1.34,  1.36,  1.38,  1.4 , 1.42]:
#        resultFiles.remove(f)

data = []
print resultFiles 

data = pyalps.loadMeasurements(resultFiles, args.y)

data = pyalps.flatten(data)
print data 

for d in data:
    d.props['observable'] =  args.y

res = pyalps.collectXY(data, x='V', y=args.y, foreach = [args.x])

print res 

print pyalps.plot.convertToText(res)
pyalps.propsort(res,args.x)

#scale data 
for d in res:
    L = d.props['L']
    d.y = [y*L**args.b for y in d.y]

    d.props['ylabel'] = r'$M_2L^{b}$'
    d.props['label'] = '$L=%g$' %(L)
    d.props['line'] = '-o'

plt.title('$b=%g$'%(args.b))

pyalps.plot.plot(res)
#plt.xlim([0,0.18])
#plt.ylim([0,0.4])

plt.legend(loc='upper left')


if args.copydata:
    for resultFile in resultFiles:
        cmd = ['cp', resultFile, '../data/']
        subprocess.check_call(cmd)

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
