'''
yL^b = f[(x-xc)*L^a]
'''

import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import os, sys , re 
import subprocess 
import socket
import argparse
from fss import FSS 
from autoScale import autoScale , myValue
from numpy import array , linspace , sqrt , argsort 

parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheaders", nargs='+', default="params", help="fileheaders")

parser.add_argument("-x", default="V", help="variable")
parser.add_argument("-y", default="M2", help="observable")
parser.add_argument("-copydata", action='store_true',  help="copy data")


group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

group = parser.add_argument_group()
group.add_argument("-Scaling", action='store_true',  help="whether to do the scaling")
group.add_argument("-FSS", action='store_true',  help="fit with fss")
group.add_argument("-Autoscale", action='store_true',  help="fit with autoscale")

group.add_argument("-xc", type = float,  default=0., help="xc")
group.add_argument("-a", type = float,  default=0., help="a")
group.add_argument("-b", type = float,  default=0., help="b")
group.add_argument("-c", type = float,  default=0., help="c")
group.add_argument("-d", type = float,  default=0., help="d")
group.add_argument("-w", type = float,  default=0., help="w")
group.add_argument("-phi", type = float,  default=0., help="phi")

#the range overwhich we perform minimization, do not mixup with xcmin and xcmax
group.add_argument("-xmin", type = float, default = -1.,  help="xmin")
group.add_argument("-xmax", type = float, default = 1.,  help="xmax")

#fix things in autoscale 
group.add_argument("-Fixxc", action='store_true',  help="fix")
group.add_argument("-Fixa",  action='store_true',  help="fix")
group.add_argument("-Fixb",   action='store_true', help="fix")
group.add_argument("-Fixc",   action='store_true', help="fix")
group.add_argument("-Fixw",   action='store_true', help="fix")


args = parser.parse_args()

resultFiles = []
for fileheader in args.fileheaders:
    resultFiles += pyalps.getResultFiles(prefix=fileheader)
print resultFiles 


#filter resultFilies
for f in list(resultFiles):
    L = int(re.search('L([0-9]*)W',f).group(1)) 
#    W = int(re.search('W([0-9]*)N',f).group(1)) 
    V= float(re.search('V([0-9]*\.?[0-9]*)ITIMEMAX',f).group(1)) 
    
    if V< 1.3 or V>1.4 or L in [3]:
        resultFiles.remove(f)



resultFiles = list(set(resultFiles))
print resultFiles 

data = pyalps.loadMeasurements(resultFiles, args.y)
data = pyalps.flatten(data)
print data 

if args.y=='nncorr':
    res = []
    for d in data:
        L = int(d.props['L'])

        r = pyalps.DataSet()
        #r.y = [d.y[-1]]
        r.y = [abs(d.y[L])]
        r.props = d.props 
        r.props['observable'] = 'farthestnncorr'
        res.append(r)

    res = pyalps.collectXY(res, x=args.x, y='farthestnncorr', foreach = ['L'])
else:
    res = pyalps.collectXY(data, x=args.x, y=args.y, foreach = ['L'])

print res

###############  FSS ############################
if args.Scaling and args.FSS:
    fss = FSS(res, xc=args.xc, a=args.a, b=args.b, c = args.c, d=args.d, w=args.w, phi = args.phi, order = 4
                 , xmin=args.xmin, xmax = args.xmax 
                 , xcmin=1.32, xcmax = 1.37
                 , amin=1.0, amax = 1.4
                 , bmin=1, bmax = 2.
                 , cmin=0., cmax = 0.
                 , dmin=0., dmax = 0.
                 , wmin=0., wmax= 2.
                 , phimin = 0., phimax = 0.
             )
    parms, chi2 = fss()

    args.xc = parms[fss.order+1]
    args.a = parms[fss.order+2]
    args.b = parms[fss.order+3]
    args.c = parms[fss.order+4]
    args.d = parms[fss.order+5]
    args.w= parms[fss.order+6] 
    args.phi = parms[fss.order+7] 

    print "xc, a, b, c, d, w, phi, chi2:", args.xc, args.a, args.b, args.c, args.d, args.w, args.phi, chi2 

############### autoScale ############################
if args.Scaling and args.Autoscale:
    data = {}  
    for r in res:
        L = int(r.props['L'])
        tmp = []
        for x, y in zip(r.x, r.y):
            tmp.append(myValue( L, x, y.mean, y.error))
        data[L] = tmp 
    print data 

    xc, a, b, c, w, S = autoScale(data, args.xc, args.a, args.b, args.c, args.w, args.xmin, args.xmax, xco=(not args.Fixxc), ao=(not args.Fixa), bo=(not args.Fixb),  co=(not args.Fixc), wo=(not args.Fixw), repFit=1, getError=1)

    print "xc, a, b, c, w, S:", xc, a, b, c, w, S 
    args.xc = xc
    args.a = a
    args.b = b
    args.c = c
    args.w = w
########################################################


if (args.Scaling):
    for d in res:
        L = d.props['L']
        d.y = [y*L**args.b/(1.+ args.c* L **(-args.w) ) for y in d.y]
        d.x = (d.x- args.xc) * L ** args.a + args.d * L **(-args.phi)

        d.props['xlabel'] = r'$(x-x_c)L^{a}+dL^{-\phi}$'
        d.props['ylabel'] = r'$yL^{b}/(1+cL^{-\omega})$'
        d.props['label'] = '$L= %s$' % int(d.props['L'])

#print pyalps.plot.convertToText(res)
pyalps.propsort(res,'L')
pyalps.plot.plot(res)
plt.legend(loc='best')

if (args.Scaling):
    plt.title('$x_c=%g,a=%g,b=%g,c=%g,d=%g,\omega=%g,\phi=%g$'%(args.xc,args.a, args.b, args.c, args.d, args.w, args.phi))


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
