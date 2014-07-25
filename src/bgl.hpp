#ifndef BGL_HPP
#define BGL_HPP

#include <iostream>
#include <algorithm> 
#include <vector>
#include <iterator>
#include <alps/lattice.h>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>

template<typename Graph, typename DistanceMap>
class my_visitor : public boost::default_bfs_visitor{

public:
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

    my_visitor(DistanceMap& dmap, const unsigned int N):dmap_(dmap), d_(N){ }

    void tree_edge(const edge_descriptor &e, const Graph &g) const {
      //std::cout << "Tree edge: " << e << std::endl;
      d_[target(e, g)] = d_[source(e, g)] + 1;
      dmap_[d_[target(e, g)]].push_back(target(e, g)); 
    }

private:
    DistanceMap& dmap_; 
    mutable std::vector<unsigned int> d_; //distances to origin  
};


std::vector<DistanceMap> get_distmap(const alps::graph_helper<>& lattice){

  typedef alps::graph_helper<>::graph_type Graph; 
     
  DistanceMap dmap; 
  my_visitor<Graph, DistanceMap> vis(dmap, lattice.num_sites());

  std::vector<DistanceMap> distmap; //each element stores a map for a site 
  for (site_type origin = 0; origin < lattice.num_sites(); ++origin){

      dmap.clear(); 
      dmap[0].push_back(origin); 
      breadth_first_search(lattice.graph(), origin, boost::visitor(vis));

      distmap.push_back(dmap);  
 }

  return distmap;
} 

Eigen::MatrixXi get_disttable(const std::vector<DistanceMap>& distmap, const unsigned int Ns){

    Eigen::MatrixXi disttable(Ns,Ns); 

    for (unsigned int si = 0; si < Ns; ++si){
      for (DistanceMap::const_iterator it = distmap[si].begin(); it != distmap[si].end(); ++it) {
        int dist = it->first; 
        for (unsigned int j=0; j<  it->second.size() ; ++j){
            site_type sj = it->second[j];
            disttable(si, sj) = dist; 
       }
      }
    }
    return disttable; 
}


std::vector<std::vector<unsigned int> >  get_neighbors(const std::vector<DistanceMap>& distmap, const unsigned int Ns, const unsigned int l)
{

  std::vector<std::vector<unsigned int> >  neighbors(Ns);  //for each site the list contains its neighbors within distances l 
  for (unsigned int s = 0; s < Ns; ++s){
   DistanceMap dmap = distmap[s]; 
   for (unsigned int i=1; i<= l; ++i) { // insert worm of the distance 1~l, no identical sites 
       std::copy (dmap[i].begin(), dmap[i].end(), std::back_inserter(neighbors[s])); 
   }
  }

  return neighbors; 
}

std::vector<unsigned int> get_shellsize(const std::vector<DistanceMap>& distmap){

    std::vector<unsigned int> shellsize; 

    for (DistanceMap::const_iterator it = distmap[0].begin(); it != distmap[0].end(); ++it) {
        shellsize.push_back(it->second.size()); 
    }
    return shellsize; 
}


std::vector<unsigned int> get_neighborshellsize(const std::vector<unsigned int>& shellsize){
    std::vector<unsigned int> neighborshellsize(shellsize.size()); 

    neighborshellsize[0] = shellsize[1]; 
    for (int dist = 1; dist < shellsize.size()-1; ++dist){
        neighborshellsize[dist] = shellsize[dist-1] + shellsize[dist+1]; 
    }
    neighborshellsize[shellsize.size()-1] = shellsize[shellsize.size()-2]; 
    
    return neighborshellsize; 
}

#endif
