#ifndef NEIGHBORS_HPP
#define NEIGHBORS_HPP

#include <vector>
#include <algorithm>

std::vector<std::vector<unsigned int> >  get_neighbors(const std::vector<DistanceMap>& distmap, const unsigned int Ns, const unsigned int l)
{

  std::vector<std::vector<unsigned int> >  neighbors(Ns);  //for each site the list contains its neighbors within distances l 
  for (unsigned int s = 0; s < Ns; ++s){
   const DistanceMap& dmap = distmap[si]; 
   for (unsigned int i=0; i<= std::min(l, dmap.size()); ++i) {
       std::copy (dmap[i].begin(), dmap[i].end(), std::back_inserter(neighbors[s])); 
   }
  }

  return neighbors; 
}
#endif
