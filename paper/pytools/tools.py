
import os, glob
import numpy as np
from itertools import chain

from dataset import DataSet
from hist import flatten

def recursiveGlob(dirname,pattern):
    ret = glob.glob(os.path.join(dirname, pattern))
    for d in os.listdir(dirname):
        d = os.path.join(dirname, d)
        if os.path.isdir(d):
            ret += recursiveGlob(d, pattern)
    return ret


def dict_intersect(dicts):
    """ computes the intersection of a list of dicts
    
        this function takes a list of dicts as input and returns a dict containing all those key-value pairs that appear with identical values in all dicts 
    """
    sets = [set(q.keys()) for q in dicts]
    intersection = sets[0]
    for iset in sets:
        intersection &= iset
    ret = {}
    for key in intersection:
        take = True
        val0 = dicts[0][key]
        for idict in dicts:
            try:
                if val0 != idict[key]:
                    take = False
            except:
                if np.all(val0 != idict[key]):
                    take = False
        if take:
            ret[key] = dicts[0][key]
    return ret


def convert_to_text(data,title=None,xaxis=None,yaxis=None):
    output = ''
    if  title!=None:
        output += title + '\n'
    
    if xaxis != None:
        output += '# X'
        if 'label' in xaxis:
            output += ': ' + xaxis['label']
        if 'min' in xaxis and 'max' in xaxis:
            output += ': ' + str(xaxis['min']) + ' to ' + str(xaxis['max'])
        output+='\n'

    if yaxis != None:
        output += '# Y'
        if 'label' in yaxis:
            output += ': ' + yaxis['label']
        if 'min' in yaxis and 'max' in yaxis:
            output += ': ' + str(yaxis['min']) + ' to ' + str(yaxis['max'])
        output+='\n\n'
    
    mydata = data
    if type(data) != list:
        mydata = [data]
    for q in flatten(mydata):
        if 'label' in q.props and q.props['label'] != 'none':
            output += '# ' + q.props['label']
        elif 'filename' in q.props:
            output += '# ' + q.props['filename']
        output += '\n'
        if 'xlabel' in q.props:
            output += '# X: ' + q.props['xlabel'] + '\n'
        if 'ylabel' in q.props:
            output += '# Y: ' + q.props['ylabel'] + '\n'

        for i in range(len(q.x)):
            output += str(q.x[i]) + '\t' + str(q.y[i]) + '\n'
        output+='\n\n'                
    return output


def collectXY(sets,x,y,foreach=[]):
      """ collects specified data from a list of DataSet objects
         
          this function is used to collect data from a list of DataSet objects, to prepare plots or evaluation. The parameters are:
    
            sets:    the list of datasets
            x:       the name of the property or measurement to be used as x-value of the collected results 
            y:       the name of the property or measurement to be used as y-value of the collected results 
            foreach: an optional list of properties used for grouping the results. A separate DataSet object is created for each unique set of values of the specified parameers.
            
          The function returns a list of DataSet objects.
      """
      foreach_sets = {}
      for iset in flatten(sets):
          if iset.props['observable'] != y:
              continue
          
          fe_par_set = tuple((iset.props[m] for m in foreach))
          
          if fe_par_set in foreach_sets:
              foreach_sets[fe_par_set].append(iset)
          else:
              foreach_sets[fe_par_set] = [iset]
      
      for k,v in foreach_sets.items():
          common_props = dict_intersect([q.props for q in v])
          res = DataSet()
          res.props = common_props
          for im in range(0,len(foreach)):
              m = foreach[im]
              res.props[m] = k[im]
          res.props['xlabel'] = x
          res.props['ylabel'] = y
          
          for data in v:
              if len(data.y)>1:
                  res.props['line'] = '.'
              xvalue = np.array([data.props[x] for i in range(len(data.y))])
              if len(res.x) > 0 and len(res.y) > 0:
                  res.x = np.concatenate((res.x, xvalue ))
                  res.y = np.concatenate((res.y, data.y))
              else:
                  res.x = xvalue
                  res.y = data.y
          
          order = np.argsort(res.x, kind = 'mergesort')
          res.x = res.x[order]
          res.y = res.y[order]
          res.props['label'] = ''
          for im in range(0,len(foreach)):
              res.props['label'] += '%s = %s ' % (foreach[im], k[im])
          
          foreach_sets[k] = res
      return foreach_sets.values()


def save_with_prov(fname, fig, obs):
    """save matplotlib figure together with provenance information"""
    dirname = 'fig_'+fname
    basename = dirname + '/' + dirname
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    fig.savefig(basename+'.pdf')
    fp = open(basename+'.txt', 'w')
    fp.write(convert_to_text(obs))
    fp.close()


def propsort(data,pn):
    '''sort datasets in data using the property named pn as key'''
    data.sort(cmp=lambda x,y:cmp(x.props[pn],y.props[pn]))
