#include "interaction_expansion.hpp"

void InteractionExpansion::measure_M2()
{
   site_type si = randomint(n_site); // we only do it for fixed si becaues of translational invarance 
    
   //itime_type itau = randomint(itime_max);
   //itime_type itau = iblock*blocksize + randomint(blocksize);// a random time inside this block 

   //Mat gtau = gf.G(itau, tlist, vlist);
   //Mat gtau = gf.wrap(itau, tlist, vlist); 

   Eigen::VectorXd denmat = gf.denmat(si); //g_{ij} for fixed i 

   //g_ji = delta_ij - eta_i eta_j * g_ij 

   //double parity = lattice.parity(si) * lattice.parity(sj);  
   //double delta = (si==sj? 1.0: 0.0); 
   //double M2 =  parity*((1. - gii)  * (1. - gjj) + (delta-gji) * gij); 
    
   double M2 =denmat.adjoint()* denmat; // g_ij * g_ji
   measurements["M2"] << M2/n_site;   
    
   double IntE = 0.0; 
   for (unsigned j=0 ; j< lattice.num_neighbors(si); ++j){
        site_type sj = lattice.neighbor(si, j);
        IntE += -denmat(sj)* denmat(sj) ; 
   }
   measurements["IntE"] << 0.5*V * IntE;// interaction energy per site : -IntE * n_site * Theta = PertOrder 

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
