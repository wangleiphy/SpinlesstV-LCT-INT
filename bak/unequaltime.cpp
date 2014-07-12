#include "interaction_expansion.hpp"
#include <vector>

void InteractionExpansion::measure_gf(){

 unsigned int Msize = M.matrix().rows();
 std::vector<double> gf(n_taumeasure); 

 for (unsigned it =0; it< n_taumeasure; ++it) {

    site_t s = (unsigned int)(random()*n_site); 
    itime_t delta_t = (double)it*beta/n_taumeasure; 

    gf[it] = green0_spline(delta_t, s, s); 
    
    if(Msize>0){
   
       Eigen::VectorXd Q(Msize);  
       Eigen::RowVectorXd R(Msize);
    
        for(unsigned int j=0; j< Msize; ++j){
             R(j) = green0_spline(delta_t-M.creators()[j].t(), s, M.creators()[j].s());
             Q(j) = green0_spline(M.creators()[j].t()-0., M.creators()[j].s(), s);
       }
    
      gf[it] -=  R * M.matrix() * Q; 
    }
  }

  measurements["gf"] << gf; 
}
 

void InteractionExpansion::measure_ntaun(){

 unsigned int Msize = M.matrix().rows();
 std::vector<double> ntaun(n_taumeasure); 

 for (unsigned it =0; it< n_taumeasure; ++it) {

    site_t s = (unsigned int)(random()*n_site); 
    itime_t delta_t = (double)it*beta/n_taumeasure; 

    std::vector<double> taus;  
    taus.push_back(delta_t);
    taus.push_back(0.);

    Eigen::Matrix2d S = Eigen::Matrix2d::Zero(); // diagonal term is always zero ;  

    //compute the Green's functions with interpolation
    S(0,1) = green0_spline(delta_t, s, s); 
    S(1,0) = -S(0,1);
    
    if(Msize>0){

        Eigen::MatrixX2d Q(Msize,2);  // M*2
        Eigen::Matrix2Xd R(2,Msize);  // 2*M   

        for(unsigned int i=0; i<2; ++i){
         for(unsigned int j=0; j< Msize; ++j){
                Q(j,i) = green0_spline(M.creators()[j].t()-taus[i], M.creators()[j].s(), s);
                R(i,j) = -lattice.parity(s) * M.creators()[j].parity()* Q(j,i);//anti-symmetrization 
          }
        }
    
      S -=  R * M.matrix() * Q; 
    }

    ntaun[it] += S.determinant(); 
  }

  measurements["ntaun"] << ntaun; 
}
