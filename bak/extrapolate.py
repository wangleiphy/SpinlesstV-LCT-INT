from pytools.floatwitherror import FloatWithError as fwe
from scipy.optimize import curve_fit 
from numpy import array, linspace , sqrt 
import matplotlib.pyplot as plt 

def extrapolate(res, nextrapolate):

    for d in res:
        # do curve_fit 
        ymean = array([y.mean for y in d.y]) # take only the mean 
        yerror = array([y.error for y in d.y]) 

        def func(x, *p):
             #return p[0]*x + p[1]
             return p[0]*x**2 + p[1]*x + p[2]
             #return p[0]*x**3 + p[1]*x**2 + p[2]*x + p[3]
        
        #try:
        popt, pcov = curve_fit(func, d.x[-nextrapolate:], ymean[-nextrapolate:], sigma = yerror[-nextrapolate], p0=array((0.,0.,0.)))
        
        xlist = linspace(0,0.1,100)
        plt.plot(xlist, func(xlist, *tuple(popt)), '--', color = d.props['color'])
        #except:
        #    pass 

        #plot intersection at 1/L = 0
        M2 = fwe()
        numbins = len(d.y[0].jackknife)
        for k in range(numbins):
            y = [yi.jackknife[k]  for yi in d.y]
            popt, pcov = curve_fit(func, d.x[-nextrapolate:], y[-nextrapolate:], p0=array((0.,0.,0.)))
            #print k, popt
            M2.jackbins.append(popt[-1])
        M2.jacknife_eval()

        if M2.mean > 0.0:
            print d.props["V"], M2  #sqrt(M2)*2.

        plt.errorbar(0, M2.mean, M2.error, linewidth=6, alpha=0.5, color =d.props['color'])
