#include "interaction_expansion.hpp"
#include <iostream>

void InteractionExpansion::print(std::ostream &os) const{
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"***                                      CTQMC for spinless fermions                                    ***"<<std::endl;
  os<<"***                                      Lei Wang, ETH Zurich, 2013                                     ***"<<std::endl;
  os<<"***                                      lewang@phys.ethz.ch                                            ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"n_site: "<<n_site << ",\tn_bond: "<< n_bond << std::endl; 
  os<<"n_tau: "<<n_tau << ",\tmax order: "<<max_order << std::endl; 
  os<<"mc steps: "<<mc_steps << ",\ttherm steps: "<<therm_steps << std::endl;
  os<<"recalc period: "<<recalc_period<<",\tmeasurement period: "<< measurement_period << std::endl; 
  os<<"T: "<<temperature<<",\tV: "<<V<<std::endl;
  
 {
  os<<"probs: ";
  std::ostream_iterator<double> out_it (os,", ");
  copy(probs.begin(), probs.end(), out_it );
 }
  os<<"\teta2, eta4: "<< eta2 << " "<< eta4 <<std::endl;

  os<<"n_max: "<< n_max <<std::endl;
  os<<"# of worm neighbors: "<< neighbors[0].size()  <<std::endl;

 { 
  std::ostream_iterator<unsigned int> out_it (os,", ");
  os<<"shellsize: ";
  std::copy (shellsize.begin(), shellsize.end(), out_it);
  os<<std::endl; 
  os<<"neighborshellsize: ";
  std::copy (neighborshellsize.begin(), neighborshellsize.end(), out_it);
  os<<std::endl; 
 }

 {
  unsigned origin = 1; 
  DistanceMap dmap = distmap[origin] ; 
  std::ostream_iterator<site_descriptor> out_it (os," ");
  for (DistanceMap::iterator it = dmap.begin(); it != dmap.end(); ++it) {
    os << origin <<  " has " << it->second.size()  << " neighbors in " << it->first << " steps: "; 
    std::copy (it->second.begin(), it->second.end(), out_it);
    os << std::endl; 
  }
 }

}


void InteractionExpansion::update_params(alps::params &parms) const{
    //std::cout << "before: " << parms.size() << std::endl; 

    parms["Nsite"] = n_site; 
    parms["Nbond"] = n_bond; 
    parms["Ncell"] = n_cell; 

    //std::cout << "after: " << parms.size() << std::endl; 
}
