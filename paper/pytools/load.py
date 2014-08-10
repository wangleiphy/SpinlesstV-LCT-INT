# Data loading tools.
# by Troyer Group, Institute for Theoretical Physics, ETH Zurich (C) 2012 - 2013
# Extended from pyalps library in ALPS 2.

import copy, os, traceback
import numpy as np

import hdf5
from dataset import DataSet

from floatwitherror import FloatWithError as fwe

# jackknife routines
import pyjack


def log(m):
    """ print a log message either to the console or through the VisTrails logger """
    print m

def parse_label(label):
    if '--' in label:
      vals = label.rsplit('--')
      ret = ()
      for val in vals:
          ret = ret + (eval(val),)
      return ret
    else:
      return eval(str(label))
 
 
def parse_labels(labels):
    if type(labels)==int: 
      return np.array([labels])
    larr=[]
    allsame = True
    first = None
    for x in labels:
      v = parse_label(x)
      larr.append(v)
      if '--' in x:      
        if first==None:
          first = v[0]
        else:
          if first != v[0]:
            allsame = False
      else:
        allsame = False
    if allsame:
      larr = [x[1] for x in larr]
    return np.array(larr)

class Hdf5Missing(Exception):
    def __init__(self,what):
        self.what = what
    
    def __str__(self):
        return 'Failed to find ' + self.what

class Hdf5Loader:
    """The Hdf5Loader class loads simulation parameters and observables from hdf5-files and returns them as hierarchical datasets"""
    def GetFileNames(self, flist):
        files = []
        for f in flist:
          if f[-4:]=='.xml':
            f = f[:-3]+'h5'
          else:
            if f[-3:]!='.h5':
              f += '.h5'
          if os.path.exists(f):
            files.append(f)
          else:
            log( "FILE "+ f+ "DOES NOT EXIST!")
        return files
        
    def ReadParameters(self,proppath):
        dict = {'filename' : self.h5fname}
        for m in self.h5f[proppath]:
                try:
                    dict[m] = self.h5f[proppath+'/'+m]
                    try:
                        dict[m] = float(dict[m])
                    except:
                        dict[m] = map(float,dict[m])
                except ValueError:
                    pass
        return dict 
        
    def GetObservableList(self,respath):
        if isinstance(self.h5f[respath], hdf5.nodeList):
            olist = self.h5f[respath]
        else:
            olist = []
        return olist
                
    # Pre: file is a hdf5 file descriptor
    # Post: returns DataSet with all parameters set
    def ReadMeasurementFromFile(self,flist,proppath='/parameters',respath='/simulation/results',measurements=None,verbose=False,jackbins=False,timeseries=False):
        fs = self.GetFileNames(flist)
        sets = []
        for f in fs:
            try:
                fileset = []
                with hdf5.archive(f) as ar:
                    self.h5f = ar
                    self.h5fname = f
                    if verbose: log("Loading from file" + f)
                    list_ = self.GetObservableList(respath)
                    params = self.ReadParameters(proppath)
                    obslist = []
                    if measurements == None:
                        obslist = list_
                    else:
                        obslist = [obs for obs in measurements if obs in list_]
                    for m in obslist:
                        if verbose: log( "Loading " + m)
                        size=0
                        xmin=0
                        xstep=1
                        if "error" in self.h5f[respath+'/'+m+'/mean']: 
                            if not isinstance(self.h5f[respath+'/'+m+'/mean/value'], np.ndarray) or not len(self.h5f[respath+'/'+m+'/mean/value'].shape):
                                mean = self.h5f[respath+'/'+m+'/mean/value']
                                error = self.h5f[respath+'/'+m+'/mean/error']
                                if jackbins and "jacknife" in self.h5f[respath+'/'+m]:
                                    jb = self.h5f[respath+'/'+m+'/jacknife/data']
                                    obs = fwe(mean, error, jackbins=jb)
                                elif jackbins and "timeseries" in self.h5f[respath+'/'+m]:
                                    ts = self.h5f[respath+'/'+m+'/timeseries/data']
                                    count = self.h5f[respath+'/'+m+'/count']
                                    binsize = int(count/len(ts))
                                    ts = ts/binsize
                                    jb = pyjack.prepare_jacknife(ts)
                                    obs = fwe(mean, error, jackbins=jb)
                                else:
                                    obs = fwe(mean, error)
                                if timeseries and "timeseries" in self.h5f[respath+'/'+m]:
                                    obs.timeseries = self.h5f[respath+'/'+m+'/timeseries/data']
                                obs=np.array([obs])
                                size=1
                            else:
                                raise Exception('Load of vector measurements not implemented.')
                        else:
                            if not isinstance(self.h5f[respath+'/'+m+'/mean/value'], np.ndarray) or not len(self.h5f[respath+'/'+m+'/mean/value'].shape):
                                obs = self.h5f[respath+'/'+m+'/mean/value']
                                obs=np.array([obs])
                                size=1
                            else:
                                obs = self.h5f[respath+'/'+m+'/mean/value']
                                size=len(obs)
                        try:
                            d = DataSet()
                            d.y = obs
                            d.x = np.arange(xmin,xmin+xstep*size,xstep)
                            d.props['hdf5_path'] = respath +"/"+ m
                            d.props['observable'] = m
                            d.props.update(params)
                            fileset.append(d)
                        except AttributeError:
                            log( "Could not create DataSet")
                    sets.append(fileset)
            except Exception, e:
                log(e)
                log(traceback.format_exc())
        return sets


def loadMeasurements(files,what=None,verbose=False,jackbins=False,timeseries=False):
    """ loads ALPS measurements from ALPS HDF5 result files
    
        this function loads results of ALPS simulations ALPS HDF5 result files
        
        The parameters are:
        
        files: a list of ALPS result files which can be either XML or HDF5 files. XML file names will be changed to the corresponding HDF5 names.
        what: an optional argument that is either a string or list of strings, specifying the names of the observables which should be loaded
        verbose: an optional boolean argument that if set to True causes more output to be printed as the data is loaded
        
        The function returns a list of list of DataSet objects. 
        The elements of the outer list each correspond to the file names specified as input.
        The elements of the inner list are each for a different observable.
        The y-values of the DataSet objects are the measurements and the x-values optionally the labels (indices) of array-valued measurements
    """
    ll = Hdf5Loader()
    if isinstance(what,str):
      what = [what]
    return ll.ReadMeasurementFromFile(files,measurements=what,verbose=verbose,jackbins=jackbins,timeseries=timeseries)

