#include "interaction_expansion.hpp"

void InteractionExpansion::measure_M2()
{

   site_type si = randomint(n_site); 
   site_type sj = randomint(n_site); 

   double gij = gf.U().row(si) * gf.gtau() * gf.Udag().col(sj); // rotate it to real space 

   double parity = lattice.parity(si) * lattice.parity(sj);  
   double gji = -parity * gij;  
   double delta = si==sj? 1.0: 0.0; 

   /*CDW structure factor*/
   double M2= parity*(delta-gji) * gij; 

   measurements["M2"] << M2; 
}
