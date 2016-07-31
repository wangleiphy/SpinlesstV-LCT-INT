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
from pytools.floatwitherror import FloatWithError as fwe
import re 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

parser.add_argument("-x", default="V", help="variable")
parser.add_argument("-y", default="M2", help="observable")
parser.add_argument("-copydata", action='store_true',  help="copy data")
parser.add_argument("-logscale", action='store_true',  help="logscale")

parser.add_argument("-extrapolate", action='store_true',  help="do curve fitting")
parser.add_argument("-nextrapolate", type=int, default= 4 ,  help="number of points used in extrapolation")


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

#filter resultFilies
#for f in list(resultFiles):
#    L = int(re.search('L([0-9]*)W',f).group(1)) 
#    V= float(re.search('V([0-9]*\.?[0-9]*)ITIMEMAX',f).group(1)) 
    
#    if V< 1.3 or V>1.4 or L in [3]:
#        resultFiles.remove(f)

#    if L in [15]:
#        resultFiles.remove(f)
    #if L in [3]:
#        resultFiles.remove(f)




data = []
print resultFiles 

if args.x == 'V' or args.x =='L':
    data = pyalps.loadMeasurements(resultFiles, args.y)
else:
    data = pyalps.loadMeasurements(resultFiles, [args.x, args.y])
#data = pyalps.flatten(data)
#print data 

#for d in data:
#    d.props['observable'] =  args.y

if args.x == 'V':
    res = pyalps.collectXY(data, x='V', y=args.y, foreach = ['L'])
    pyalps.propsort(res,'L')
elif args.x == 'L':
    res = pyalps.collectXY(data, x='L', y=args.y, foreach = ['V'])
    pyalps.propsort(res,'V')
else:
    res = pyalps.ResultsToXY(data, args.x, args.y, foreach = ['V']) 

print pyalps.plot.convertToText(res)

#scale data 
if args.x =='V':
    for d in res:
        L = d.props['L']
        d.x = (d.x- args.xc) * L ** args.a 
        d.y = array([y*L**args.b for y in d.y])

        d.props['xlabel'] = r'$(V-V_c)L^{a}$'
        d.props['ylabel'] = r'$M_2L^{b}$'
        d.props['label'] = '$L=%g$' %(L)
        d.props['line'] = '-o'

        if args.y == 'PertOrder':
            d.y /= (2.*L**2)

elif args.x == 'L':
    for d in res:
        d.x = 1./d.x 
        d.props['xlabel'] = '$1/L$'
        d.props['line'] = '-o'
else:
    pass 

if args.extrapolate:
    for d in res:
        # do curve_fit 
        from scipy.optimize import curve_fit 
        ymean = array([y.mean for y in d.y]) # take only the mean 
        yerror = array([y.error for y in d.y]) 

        def func(x, *p):
             return p[0]*x + p[1]
             #return p[0]*x**2 + p[1]*x + p[2]
             #return p[0]*x**3 + p[1]*x**2 + p[2]*x + p[3]
        
        #try:
        popt, pcov = curve_fit(func, d.x[-args.nextrapolate:], ymean[-args.nextrapolate:], sigma = yerror[-args.nextrapolate], p0=array((0.,0.)))
        
        print d.props['V'], popt
        
        xlist = linspace(0,0.1,100)
        plt.plot(xlist, func(xlist, *tuple(popt)), '--')
        #except:
        #    pass 


        #plot intersection at 1/L = 0
        M2 = fwe()
        numbins = len(d.y[0].jackknife)
        for k in range(numbins):
            y = [yi.jackknife[k]  for yi in d.y]
            popt, pcov = curve_fit(func, d.x[-args.nextrapolate:], y[-args.nextrapolate:], p0=array((0.,0.)))
            #print k, popt
            M2.jackbins.append(popt[1])
        M2.jacknife_eval()

        plt.errorbar(0, M2.mean, M2.error, linewidth=6, alpha=0.5)


plt.title('$V_c, a, b = %g, %g, %g$'%(args.xc, args.a, args.b))

pyalps.plot.plot(res)
#plt.xlim([1.28,1.42])
#plt.ylim([0,0.4])

plt.legend(loc='best')

if args.logscale:
    plt.gca().set_yscale('log')

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
