#include "interaction_expansion.hpp"

void InteractionExpansion::add()
{
  //std::cout << "##add##" << std::endl; 

  if(tlist.size()+1 > max_order) 
    return; 
    
  tlist_type::const_iterator lower, upper; 
  lower = std::lower_bound (tlist.begin(), tlist.end(), iblock*blocksize); 
  upper = std::upper_bound (tlist.begin(), tlist.end(), (iblock+1)*blocksize, std::less_equal<itime_type>());  //equal is exclude

  unsigned num_vertices = std::distance(lower, upper);  //number of vertices in this block

  itime_type itau = iblock*blocksize + randomint(blocksize);// a random time inside this block 

  if (tlist.find(itau) != tlist.end()) // there is already a vertex at itau, do nothing: we can not have two vertex at the same tau 
      return; 

  std::vector<site_type> sites; 

  alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
  sites.push_back(lattice.source(b));
  sites.push_back(lattice.target(b));
    
  // true means compute_only_weight
  double metropolis_weight = -0.25*(beta/nblock)*V*n_bond/(num_vertices+1)*add_impl(itau, sites, true) * (Remove/Add) ;


  if(fabs(metropolis_weight) > random()){

    measurements["Add"] << 1.;
    add_impl(itau, sites, false);


    sign*=metropolis_weight<0.?-1.:1.;
  }else{

    measurements["Add"] << 0.;

  }
}


void InteractionExpansion::remove()
{
    //std::cout << "##remove##" << std::endl; 
    if(tlist.size()< 1)
      return;    

    tlist_type::const_iterator lower, upper; 
    lower = std::lower_bound (tlist.begin(), tlist.end(), iblock*blocksize); 
    upper = std::upper_bound (tlist.begin(), tlist.end(), (iblock+1)*blocksize, std::less_equal<itime_type>());  //equal is exclude

    unsigned num_vertices = std::distance(lower, upper); //number of vertices in this block

    if(num_vertices < 1){
        return; 
    }
    
    std::advance(lower, randomint(num_vertices)); //the vertex to remove 
    itime_type itau = *lower; 

    //std::cout << "itau to remove " << itau << std::endl; 

    double metropolis_weight = 4.*num_vertices/(-(beta/nblock)*V*n_bond) * remove_impl(itau, true) *  (Add/Remove);

    if(fabs(metropolis_weight) > random()){ //do the actual update

      measurements["Removal"] << 1.;
      remove_impl(itau, false);  // false means really perform, not only compute weight

      sign*=metropolis_weight<0.?-1.:1.;

    }else{

      measurements["Removal"] << 0.;

    }
}

void InteractionExpansion::shift()
{

  if(tlist.empty()) 
    return; 
  
  tlist_type::const_iterator lower, upper; 
  lower = std::lower_bound (tlist.begin(), tlist.end(), iblock*blocksize); 
  upper = std::upper_bound (tlist.begin(), tlist.end(), (iblock+1)*blocksize, std::less_equal<itime_type>());  //equal is exclude

  unsigned num_vertices = std::distance(lower, upper); //number of vertices in this block

  // no vertex in current block do nothing
  if(num_vertices == 0){
      return; 
  }

  std::advance(lower, randomint(num_vertices)); //the vertex to shift 
  itime_type itau = *lower; 
    
  //mv sj to sjprime, which is another neighbor of si (but different from sj)
  unsigned i = randomint(2); // 0 or 1 
  site_type si = vlist[itau][i]; 
  site_type sj = vlist[itau][1-i]; 

  //std::cout << "...2" << std::endl; 

  std::vector<site_type> neighbors; 
  for (unsigned j=0 ; j< lattice.num_neighbors(si); ++j){
       if (lattice.neighbor(si, j)!=sj){
            neighbors.push_back(lattice.neighbor(si, j)); 
       }
  }

  //std::cout << "sj, neighbors " <<  sj << " " << neighbors[0] << " " << neighbors[1] << std::endl; 

  std::vector<site_type> sites; 
  sites.push_back(si); 
  sites.push_back(sj);
  sites.push_back(neighbors[randomint(neighbors.size())] ); //sjprime 
    
  // true means compute_only_weight
  double metropolis_weight = shift_impl(itau, sites, true);

  //std::cout << "...4 " <<  metropolis_weight << std::endl; 

  if(fabs(metropolis_weight) > random()){

    measurements["Shift"] << 1.;
    shift_impl(itau, sites, false);

    sign*=metropolis_weight<0.?-1.:1.;
  }else{

    measurements["Shift"] << 0.;

  }
}
