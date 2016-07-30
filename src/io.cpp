#include "interaction_expansion.hpp"
#include <iostream>

void InteractionExpansion::print(std::ostream &os) const{
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"***                                      LCT-INT for spinless fermions                                  ***"<<std::endl;
  os<<"***                                      Lei Wang, ETH Zurich, 2014-2015                                ***"<<std::endl;
  os<<"***                                                <lewang@phys.ethz.ch>                                ***"<<std::endl;
  os<<"***                                                IOP/CAS, 2016                                        ***"<<std::endl;
  os<<"***                                                <wanglei@iphy.ac.cn>                                 ***"<<std::endl;
  os<<"***                                      Code used in PRB 91, 235151 (2015)                             ***"<<std::endl;
  os<<"***                                      http://dx.doi.org/10.1103/PhysRevB.91.235151                   ***"<<std::endl;
  os<<"***********************************************************************************************************"<<std::endl;
  os<<"n_site: " << n_site << ",\tn_bond: "<< n_bond << std::endl; 
  os<<"max order: "<<max_order << std::endl; 
  os<<"mc steps: "<<mc_steps << ",\ttherm steps: "<<therm_steps << std::endl;
  os<<"itime_max: " <<  itime_max << ",\tnblock: "<<nblock<<",\tsteps_per_block: "<< steps_per_block << std::endl; 
  os<<"recalc_period: "<< recalc_period<<",\tblocksize: "<< blocksize << std::endl; 
  os<<"beta: "<<beta<<",\tV: "<<V<<std::endl;

 { 
  std::ostream_iterator<unsigned> out_it (os,", ");
  os<<"shellsize: ";
  std::copy (shellsize.begin(), shellsize.end(), out_it);
  os<<std::endl; 
 }
}
