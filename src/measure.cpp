#include "interaction_expansion.hpp"

void InteractionExpansion::measure_M2()
{

   itime_type itau = randomint(itime_max);
   Mat gtau = gf.G(itau, tlist, vlist); // gf at the current time 

   /*CDW structure factor*/
   double M2 = 0.; 
   for (site_type i =0; i< n_site; ++i){  
          for (site_type j =0; j< n_site; ++j){  
           double delta = i==j? 1.0: 0.0; 
           M2 += lattice.parity(i) * lattice.parity(j) * 
                 ((1. - gtau(i,i))  * (1. - gtau(j,j)) + (delta-gtau(j,i)) * gtau(i,j)); 
          }
   }

   measurements["M2"] << M2/(n_site*n_site); 
}
