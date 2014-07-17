#include "interaction_expansion.hpp"

void InteractionExpansion::add()
{
  //std::cout << "##add##" << std::endl; 

  if(tlist.size()+1 > max_order) 
    return; 
    
  tlist_type::const_iterator lower, upper; 
  lower = std::lower_bound (tlist.begin(), tlist.end(), iblock*blocksize); 
  upper = std::upper_bound (tlist.begin(), tlist.end(), (iblock+1)*blocksize, std::less_equal<itime_type>());  //equal is exclude

  //std::cout << "0" << std::endl; 

  unsigned num_vertices = std::distance(lower, upper);  //number of vertices in this block

  itime_type itau = iblock*blocksize + randomint(blocksize);// a random time inside this block 

  if (tlist.find(itau) != tlist.end()) // we can not have two vertex at the same tau 
      return; 

  //std::cout << "1" << std::endl; 

  std::vector<site_type> sites; 

  alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
  sites.push_back(lattice.source(b));
  sites.push_back(lattice.target(b));
    
  //std::cout << "2" << std::endl; 

  // true means compute_only_weight
  double metropolis_weight = -0.25*(beta/nblock)*V*n_bond/(num_vertices+1)*add_impl(itau, sites, true);

  //std::cout << "3" << std::endl; 

  if(fabs(metropolis_weight) > random()){

    measurements["Add"] << 1.;
    add_impl(itau, sites, false);

  //std::cout << "4" << std::endl; 

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
    
    //std::cout << "iblock, bounds: " << iblock << " " << iblock*blocksize << " " << (iblock+1)*blocksize << std::endl; 
    //std::cout << "vertices number " << num_vertices << std::endl; 
    //std::copy(lower, upper, std::ostream_iterator<itime_type>(std::cout, " "));
    //std::cout << std::endl;  

    std::advance(lower, randomint(num_vertices)); //the vertex to remove 
    itime_type itau = *lower; 

    //std::cout << "itau to remove " << itau << std::endl; 

    double metropolis_weight = 4.*num_vertices/(-(beta/nblock)*V*n_bond) * remove_impl(itau, true);

    if(fabs(metropolis_weight) > random()){ //do the actual update

      measurements["Removal"] << 1.;
      remove_impl(itau, false);  // false means really perform, not only compute weight

      //std::cout << "remove " << n << " vertices." << "creators: ";  
      //for (unsigned int  i=0; i< M.creators().size(); ++i) {
      //  std::cout << M.creators()[i].s()<< "("<< M.creators()[i].t() << ")"  << ","; 
      //}
      //std::cout << std::endl; 
      sign*=metropolis_weight<0.?-1.:1.;

    }else{

      measurements["Removal"] << 0.;

    }
}
