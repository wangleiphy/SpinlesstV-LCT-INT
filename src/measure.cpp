#include "interaction_expansion.hpp"

void InteractionExpansion::measure_M2()
{

   site_type si = randomint(n_site); 
   site_type sj = randomint(n_site); 

   double gij = gf.U().row(si) * gf.gtau() * gf.Udag().col(sj); // rotate it to real space 

   //g_ji = delta_ij - eta_i eta_j * g_ij 

   measurements["M2"] <<  (gij * gij); 
}
