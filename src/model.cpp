#include "interaction_expansion.hpp"

void InteractionExpansion::add()
{
  //std::cout << "##add##" << std::endl; 

  int pert_order = tlist.size(); 

  if(pert_order+1 > max_order) 
    return; 

  itime_type itau = randomint(itime_max);

  if (tlist.find(itau) != tlist.end()) // we can not have two vertex at the same tau 
      return; 

  std::vector<site_type> sites; 

  alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
  sites.push_back(lattice.source(b));
  sites.push_back(lattice.target(b));


  // true means compute_only_weight
  double metropolis_weight = -0.25*beta*V*n_bond/(pert_order+1)*add_impl(itau, sites, true);

  //if (metropolis_weight<0. && fabs(metropolis_weight) > 1E-10){
  //  std::cout << metropolis_weight << " < 0 in add" << std::endl; 

  //  std::cout << "itau, b:"  << itau << " " << b << std::endl; 
  //  std::cout << "tlist: "; 
  //  std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
  //  std::cout << std::endl; 
  //  abort();  
  //}

  if(fabs(metropolis_weight) > random()){

    std::stringstream obs_name;
    measurements["Add"] << 1.;

    add_impl(itau, sites, false);

    //std::cout << "add " << n << " vertices." << "creators: ";  
    //for (unsigned int  i=0; i< M.creators().size(); ++i) {
    //  std::cout << M.creators()[i].s()<< "("<< M.creators()[i].t() << ")"  << ","; 
    //}
    //std::cout << std::endl; 
    sign*=metropolis_weight<0.?-1.:1.;
  }else{

    measurements["Add"] << 0.;

  }
}


void InteractionExpansion::remove()
{
    //std::cout << "##remove##" << std::endl; 
    unsigned pert_order = tlist.size(); 

    if(pert_order < 1)
      return;    

    unsigned vertex_nr=randomint(pert_order);  

    double metropolis_weight = 4.*pert_order/(-beta*V*n_bond) * remove_impl(vertex_nr, true);
    //std::cout << "after remove_impl" << std::endl; 
    //if (metropolis_weight<0.){
    //  std::cout << metropolis_weight << " < 0 in remove" << std::endl; 
    // abort();  
    //}

    if(fabs(metropolis_weight) > random()){ //do the actual update

      measurements["Removal"] << 1.;
      remove_impl(vertex_nr, false);  // false means really perform, not only compute weight

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
