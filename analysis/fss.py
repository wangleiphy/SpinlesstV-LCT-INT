from scipy.optimize import fmin_l_bfgs_b
from numpy import array ,polyfit , poly1d 

'''
y *L^b = (1+cL^{-w}) f((x-xc)L^a+d*L^{-phi})
'''

class FSS:
    def __init__(self, data, xc=0.5, a=1., b=1., c= 0., d=0., w= 0., phi=0., order=4
                , xmin = -1, xmax = 1.
                , xcmin=0.5, xcmax = 0.6
                , amin=0, amax = 2. 
                , bmin=0, bmax = 2.
                , cmin=0., cmax = 1.
                , dmin=-1., dmax = 0.
                , wmin=0, wmax= 2.
                , phimin = 0., phimax = 2.
                ):
        self.data = data 
        self.order = order 
        
        #overwhich we perform optimization 
        self.xmin, self.xmax = xmin, xmax 

        #use one of the dataset to get initial parms of the polynominal 
        d = self.data[-1]
        L = d.props['L']
        x = (d.x-xc)* L**a + c * L**(-phi)
        y = [d.y[i].mean* L**b/(1.+c*L**(-w)) for i in range(len(d.y))]
        self.iniparms = list(polyfit(x, y, self.order))

        #print self.iniparms
        #p = poly1d(self.iniparms)
        #xp = linspace(-2, 2, 100)
        #plt.plot(x, y, '.', xp, p(xp), '-')
        #plt.show()

        #self.iniparms = [0.]*(1+self.order)
        self.iniparms += [xc, a, b, c, d, w, phi] # initial params of critical exponents 

        self.bounds = [(-5., 5.)] * (1+order) #boundary for the coef of poly 
        self.bounds += [(xcmin, xcmax)] # xc 
        self.bounds += [(amin, amax)]    # a 
        self.bounds += [(bmin, bmax)]   # b
        self.bounds += [(cmin, cmax)] #c 
        self.bounds += [(dmin, dmax)] #d
        self.bounds += [(wmin, wmax)]  # w
        self.bounds += [(phimin, phimax)]  #phi 

        #print len(self.iniparms), len(self.bounds)

    def g(self, coef, x):
        assert (len(coef) == self.order+1)

        res = 0.0
        for i in range(self.order+1):
            res += coef[i] * x** (self.order-i)
        return res 

    def func(self, parms):
        '''
        target function 
        '''
        assert(len(parms)==self.order +1 + 7 )

        xc = parms[self.order+1]
        a = parms[self.order+2]
        b = parms[self.order+3]
        c = parms[self.order+4]
        d = parms[self.order+5]
        w= parms[self.order+6] 
        phi = parms[self.order+7] 

        res = 0. 
        for data in self.data:
            for i in range(len(data.x)):
                L = data.props['L']
                x = (data.x[i]-xc)* L**a + d * L**(-phi)
                if x > self.xmin and x< self.xmax:
                    res += (data.y[i].mean - L**(-b)*(1.+c*L**(-w))*self.g(parms[:self.order+1], x) )**2/data.y[i].error **2
        return res 


    def __call__(self):

        res = fmin_l_bfgs_b(self.func, self.iniparms ,approx_grad=True, bounds=self.bounds)

        return res[0], res[1]


if __name__=='__main__':
    import pyalps 
    import pyalps.plot 
    from glob import glob 
    from numpy import loadtxt , linspace 
    from pyalps.alea import MCScalarData 
    import re 
    import matplotlib.pyplot as plt 
    
    #read in data 
    data = []
    for filename in glob('data/*.dat'):
        x, y, yerr = loadtxt(filename, unpack=True)

        d = pyalps.DataSet()
        d.x = x
        d.y = array([MCScalarData(mean, error) for mean, error in zip(y, yerr)])
        L = int(re.search('L([0-9]*).dat',filename).group(1)) 
        d.props['L'] = L 

        data.append(d)
    
    print data 
    #plot raw data  
    #pyalps.plot.plot(data)
    #plt.show()

    fss = FSS(data, xc=0.5927, a=0.75, b=0.104, c = 0., d=0., w= 0., phi=0., order = 4)
    parms, chi2 = fss()
     
    xc = parms[self.order+1]
    a = parms[self.order+2]
    b = parms[self.order+3]
    c = parms[self.order+4]
    d = parms[self.order+5]
    w= parms[self.order+6] 
    phi = parms[self.order+7] 

    #print parms[:fss.order]
    print "xc, a, b,c,d, w, phi, chi2:", xc, a, b, c, d, w, phi, chi2 
    #scaling according to solution 
    for dd in data:
        L = dd.props['L']
        dd.x = (dd.x-xc)* L**a + d* L**(-phi) 
        dd.y = [y * L**b/(1.+c*L**(-w))  for y in dd.y]
        #d.props['line'] = ''
    
    #the universal scaling function 
    xlist = linspace(-2, 2, 100)
    glist = [fss.g(parms[:fss.order+1], x) for x in xlist]
    plt.plot(xlist , glist, 'k-', lw =2)
    
    #scaled data 
    pyalps.plot.plot(data)
    plt.xlabel('$(x-x_c)L^a+dL^{-\phi}$')
    plt.ylabel('$yL^b/(1+cL^{-\omega})$')
    plt.xlim([-2,2])
    plt.ylim([0,1])

    plt.show()
