#include "interaction_expansion.hpp"
#include <iostream>

void InteractionExpansion::print(std::ostream &os) const{
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"***                                      CTQMC for spinless fermions                                    ***"<<std::endl;
  os<<"***                                      Lei Wang, ETH Zurich, 2013                                     ***"<<std::endl;
  os<<"***                                      lewang@phys.ethz.ch                                            ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"n_bond: "<< n_bond << std::endl; 
  os<<"max order: "<<max_order << std::endl; 
  os<<"mc steps: "<<mc_steps << ",\ttherm steps: "<<therm_steps << std::endl;
  os<<"recalc period: "<<recalc_period<<",\tmeasurement period: "<< measurement_period << std::endl; 
  os<<"T: "<<temperature<<",\tV: "<<V<<std::endl;
}

