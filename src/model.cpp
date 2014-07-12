#include "interaction_expansion.hpp"

void InteractionExpansion::add()
{
  //std::cout << "##add##" << std::endl; 

  int pert_order = M.num_vertices(); 

  if(pert_order+1 > max_order) 
    return; 

  double tau = beta*random();

  std::vector<site_t> sites; 

  alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
  sites.push_back(lattice.source(b));
  sites.push_back(lattice.target(b));

  // true means compute_only_weight
  double metropolis_weight = -0.25*beta*V*n_bond/(pert_order+1)*add_impl(tau, sites, true);

  //if (metropolis_weight<0.){
  //  std::cout << metropolis_weight << " < 0 in add" << std::endl; 
  //  abort();  
  //}

  if(fabs(metropolis_weight) > random()){

    std::stringstream obs_name;
    measurements["Add"] << 1.;

    add_impl(tau, sites, false);

    //std::cout << "add " << n << " vertices." << "creators: ";  
    //for (unsigned int  i=0; i< M.creators().size(); ++i) {
    //  std::cout << M.creators()[i].s()<< "("<< M.creators()[i].t() << ")"  << ","; 
    //}
    //std::cout << std::endl; 
 
    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 

    sign*=metropolis_weight<0.?-1.:1.;
  }else{

    measurements["Add"] << 0.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 
  }
}


void InteractionExpansion::remove()
{
    //std::cout << "##remove##" << std::endl; 
    unsigned pert_order = M.num_vertices(); 

    if(pert_order < 1)
      return;    

    unsigned vertex_nr=randomint( pert_order)

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

      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()); 

      sign*=metropolis_weight<0.?-1.:1.;

    }else{

      measurements["Removal"] << 0.;

      //do nothing
      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()); 
    }
}
