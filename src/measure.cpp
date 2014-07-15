#include "interaction_expansion.hpp"

void InteractionExpansion::measure_M2()
{
   site_type si = randomint(n_site); // we only do it for fixed si becaues of translational invarance 
    
   //itime_type itau = randomint(itime_max);
   //itime_type itau = iblock*blocksize + randomint(blocksize);// a random time inside this block 

   //Mat gtau = gf.G(itau, tlist, vlist);
   //Mat gtau = gf.wrap(itau, tlist, vlist); 

   Eigen::VectorXd denmat = gf.U() * (gf.gtau() * gf.Udag().col(si)); // rotate it to real space, N^2 operation

   //g_ji = delta_ij - eta_i eta_j * g_ij 

   //double parity = lattice.parity(si) * lattice.parity(sj);  
   //double delta = (si==sj? 1.0: 0.0); 
   //double M2 =  parity*((1. - gii)  * (1. - gjj) + (delta-gji) * gij); 
    
   double M2 =denmat.adjoint()* denmat; // g_ij * g_ji
   measurements["M2"] << M2/n_site;   
}

void InteractionExpansion::measure_vhist(){
    
    std::vector<double> vhist(nblock);  

    for (unsigned i=0; i< nblock; ++i) {

        tlist_type::const_iterator lower, upper; 
        lower = std::lower_bound (tlist.begin(), tlist.end(), i*blocksize); 
        upper = std::upper_bound (tlist.begin(), tlist.end(), (i+1)*blocksize, std::less_equal<itime_type>());  //equal is exclude
 
        vhist[i] = 1.*std::distance(lower, upper); //number of vertices in this block

    }
    measurements["Vhist"] << vhist; 
}
