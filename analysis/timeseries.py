import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
import os , sys

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument("-fileheader",default="params", help="fileheader")
parser.add_argument("-y", default="PertOrder", help="observable")

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-show", action='store_true',  help="show figure right now")
group.add_argument("-outname", default="result.pdf",  help="output pdf file")

args = parser.parse_args()

resultFiles = pyalps.getResultFiles(prefix=args.fileheader)
timeseries = pyalps.loadTimeSeries(resultFiles[0], args.y)

if len(timeseries.shape) ==2: 
    timeseries = timeseries[:,-1]

plt.figure()
plt.xlabel('t')
plt.ylabel(args.y)
plt.plot(timeseries,'-o')

plt.figure()
n, bins, patches = plt.hist(timeseries, range=(timeseries.min(), timeseries.max()), bins =20)


if args.show:
    plt.show()
else:
    plt.savefig(args.outname, dpi=300, transparent=True)
    
    #email it to me 
    pyalps.sendmail('lewang@phys.ethz.ch'    # email address of recipients 
            , message='Send from ' + os.getcwd() + ' with python ' + ' '.join([str(a) for a in sys.argv])
                   , attachment= args.outname 
                   , subject='Figure: ' + args.outname 
                   )
    os.system('rm '+args.outname)
