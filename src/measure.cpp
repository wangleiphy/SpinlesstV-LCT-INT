#include "interaction_expansion.hpp"

void InteractionExpansion::measure_M2()
{
    
   itime_type itau = randomint(itime_max);

   site_type si = randomint(n_site); 
   site_type sj = randomint(n_site); 

   //Mat gtau = gf.G(itau, tlist, vlist);

   Mat gtau = gf.wrap(itau, tlist, vlist);  

   double gij = gf.U().row(si) * gtau * gf.Udag().col(sj); // rotate it to real space 


   //g_ji = delta_ij - eta_i eta_j * g_ij 

   //double parity = lattice.parity(si) * lattice.parity(sj);  
   //double delta = (si==sj? 1.0: 0.0); 
   //double M2 =  parity*((1. - gii)  * (1. - gjj) + (delta-gji) * gij); 

   measurements["M2"] <<  gij * gij; 
}
