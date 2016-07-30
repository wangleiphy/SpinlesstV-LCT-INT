#include "interaction_expansion.hpp"

/*
void InteractionExpansion::measure_M2()
{

     Mat gtau= gf.gtau(); 

     double M2 = 0.0 ;  
      for (site_type si=0; si< n_site; ++si){
        for (site_type sj=0; sj< n_site; ++sj){
            double delta = (si==sj)? 1.0: 0.0; 
            M2 += (delta-gtau(si, sj)) * gtau(sj, si) * lattice.parity(si) * lattice.parity(sj);  
        }
      }
     measurements["M2"] << M2/(n_site*n_site);   

     double IntE = 0.0; 
     for (unsigned b =0; b< n_bond; ++b){ // loop over bonds to calculate interacion energy 
            alps::graph_helper<>::bond_descriptor bond= lattice.bond(b);
            site_type si = lattice.source(bond);
            site_type sj = lattice.target(bond); 

            IntE += -gtau(sj, si) * gtau(si, sj); 
       }
       measurements["IntE"] << V*IntE/n_site;// interaction energy per site : -IntE * n_site * Theta = PertOrder 

}
*/

void InteractionExpansion::measure_M2()
{
   site_type si = randomint(n_site); // we only do it for fixed si becaues of translational invarance 
    
   //Vec denmat = gf.denmat(si); 
   Vec denmat = gf.denmathalfTheta(si, tlist, vlist); //g_{ij}(Theta/2) for fixed si 

   //g_ji = delta_ij - eta_i eta_j * g_ij 

   //double parity = lattice.parity(si) * lattice.parity(sj);  
   //double delta = (si==sj? 1.0: 0.0); 
   //double M2 =  parity*((1. - gii)  * (1. - gjj) + (delta-gji) * gij); 
    
   double M2 =denmat.adjoint()* denmat; // g_ij * g_ji
   measurements["M2"] << M2/n_site;   
   measurements["M2k"] << (M2/n_site) * double(tlist.size());

//   double Kappa = 0.; 
//   for (unsigned sj=0 ; sj< n_site; ++sj){
//       Kappa += lattice.parity(si) * lattice.parity(sj) * denmat(sj)*denmat(sj);
//   }
//   measurements["Kappa"] << Kappa/beta/n_site; 

    
   double IntE = 0.0; 
   double KinE = 0.0; 
   for (unsigned j=0 ; j< lattice.num_neighbors(si); ++j){
        site_type sj = lattice.neighbor(si, j);
        IntE += -denmat(sj)* denmat(sj) ; 
        KinE += -K_(si, sj)* denmat(sj) ; 
   }
   measurements["IntE"] << 0.5*V * IntE;// interaction energy per site : -IntE * n_site * Theta = PertOrder 
   measurements["KinE"] << KinE; //kinetic energy per site , no 0.5 because we have <c_i^\dagger c_j > + <c_j^\dagger c_i>
   measurements["Energy"]  <<  0.5*V *IntE + KinE; // total enegy per site  
   
   measurements["KinEk"] << KinE * double(tlist.size());

  //C(R) = <(n_{i}-0.5) (n_{i+R}-0.5)> , average over site i and sites within each R shell 
  std::vector<double> nncorr(shellsize.size()); 
   for (site_type sj=0 ; sj< n_site; ++sj){
    unsigned dist = disttable(si,sj); 
    double parity = lattice.parity(si) * lattice.parity(sj); 
    nncorr[dist] +=  parity * denmat(sj) * denmat(sj) /shellsize[dist]; // normalize by number of sites in this shell 
   }
   measurements["nncorr"]  << nncorr; 

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

void InteractionExpansion::measure_M4(){

         Mat denmat = gf.denmathalfTheta(tlist, vlist); //g_{ij}(itau) 

         //measure M4 
         double M4 = 0.0; 
         site_type i = randomint(n_site); 
         for (site_type j=0; j< n_site; ++j){
             double delta_ij = (i==j)? 1.0:0.0; 
             for (site_type k=0; k< n_site; ++k){
                 double delta_ki = (k==i)? 1.0:0.0; 
                 double delta_kj = (k==j)? 1.0:0.0; 
                 for (site_type l=0; l< n_site; ++l){
         
                       Eigen::Matrix4d mat44=Eigen::Matrix4d::Zero();  
         
                       //mat44(0, 0) = denmat(i, i)-1.0; 
                       //mat44(1, 1) = denmat(j, j)-1.0; 
                       //mat44(2, 2) = denmat(k, k)-1.0; 
                       //mat44(3, 3) = denmat(l, l)-1.0; 
         
                       mat44(0, 1) = denmat(i, j); 
                       mat44(0, 2) = denmat(i, k); 
                       mat44(0, 3) = denmat(i, l); 
                       mat44(1, 2) = denmat(j, k); 
                       mat44(1, 3) = denmat(j, l); 
                       mat44(2, 3) = denmat(k, l); 
                        
                       mat44(1,0) = denmat(j, i) - delta_ij;   
                       mat44(2,0) = denmat(k, i) - delta_ki;   
                       mat44(2,1) = denmat(k, j) - delta_kj;   
                       mat44(3,0) = denmat(l, i) - ((l==i)? 1.0:0.0);   
                       mat44(3,1) = denmat(l, j) - ((l==j)? 1.0:0.0);   
                       mat44(3,2) = denmat(l, k) - ((l==k)? 1.0:0.0);   
         
                       //std::cout <<  mat44.determinant()  << std::endl; 
         
                       double parity = lattice.parity(i)*lattice.parity(j)*lattice.parity(k)*lattice.parity(l); 
                       M4 += parity * mat44.determinant(); 
                }
             }
         }
         measurements["M4"] << M4/std::pow(double(n_site), 3);  
}

