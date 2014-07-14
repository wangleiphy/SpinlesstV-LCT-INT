import pyalps.hdf5 

def report(h5_outfile, observables):
  if isinstance(observables, str):
    observables = [observables];

  results = [];
  for observable in observables:
    measurements = pyalps.hdf5.archive(h5_outfile, 'r')["/simulation/results/" + observable];

    mean  = measurements['mean']['value'];
    error = measurements['mean']['error'];
    tau   = measurements['tau']['value'];
    count = measurements['count'];
    results.append({'observable': observable,'mean': mean, 'error': error, 'tau': tau, 'count': count});

  #return results;
  for d in results:
      print d['observable']  , d['mean'] , '+/-',  d['error'] , d['tau'], d['count']
