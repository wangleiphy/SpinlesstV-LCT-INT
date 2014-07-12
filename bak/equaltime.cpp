#include "interaction_expansion.hpp"
#include <alps/multi_array.hpp>
#include <boost/foreach.hpp>

#include <set>
#include <vector>
#include <algorithm> 

void InteractionExpansion::measure_equaltime()
{

  typedef Eigen::Matrix< double , Eigen::Dynamic, Eigen::Dynamic > Mat; 
  typedef Eigen::Matrix< double , 1, Eigen::Dynamic >  RowVector; 
  typedef Eigen::Matrix< double , Eigen::Dynamic, 1 >  ColVector; 

  /*equal time green function*/
  unsigned int si= (unsigned int)(random()*n_site);  // pick up a random site as origin

  RowVector ci_cjdag = bare_green_itime.equaltime().row(si); // <c_i c^{+}_j> 
  ColVector cj_cidag = bare_green_itime.equaltime().col(si); // <c_j c^{+}_i> 

  //for (unsigned int sj=0; sj<n_site; ++sj) {             
  //    ci_cjdag(sj) = bare_green_itime.value(0, si, sj); 
  //    cj_cidag(sj) = bare_green_itime.value(0, sj, si); 
  //}
  //
  //std::cout << "cj_cidag\n" << cj_cidag << std::endl; 
  //std::cout << "ci_cjdag\n" << ci_cjdag << std::endl; 
  //std::cout << "diff\n" << ci_cjdag-cj_cidag.transpose() << std::endl; 
  //abort(); 

  if (M.matrix().rows()>0){
       double tau = beta*random();

       {//corrections to ci_cjdag 

           RowVector gf_left(M.matrix().rows()); 
           Mat gf_right(M.matrix().rows(), n_site); 
           for (unsigned int p=0; p<M.matrix().rows();++p) {
             gf_left(p) = green0_spline(tau-M.creators()[p].t(), si, M.creators()[p].s());
        
             for (unsigned int sj=0; sj<n_site; ++sj) {             
                gf_right(p, sj) = green0_spline(M.creators()[p].t()-tau, M.creators()[p].s(), sj);
             }
           }
           ci_cjdag -= (gf_left* M.matrix()) * gf_right; 
       }

       {//corrections to cj_cidag 

           Mat gf_left(n_site, M.matrix().rows()); 
           ColVector gf_right(M.matrix().rows()); 
           for (unsigned int p=0; p<M.matrix().rows();++p) {
        
              for (unsigned int sj=0; sj<n_site; ++sj) {             
                 gf_left(sj, p) = green0_spline(tau-M.creators()[p].t(), sj, M.creators()[p].s());
              }
              gf_right(p) = green0_spline(M.creators()[p].t()-tau, M.creators()[p].s(), si);
           }
           cj_cidag -= gf_left* (M.matrix() * gf_right); 
       }
  }

  //ColVector cj_cidag(n_site); // <c_j c^{+}_i> 
  //for (unsigned int sj=0; sj<n_site; ++sj) {             
  //    cj_cidag(sj) = lattice.parity(sj) * lattice.parity(si) * ci_cjdag(sj);
  //}

   //std::cout << "cj_cidag\n" << cj_cidag << std::endl; 
   //std::cout << "ci_cjdag\n" << ci_cjdag << std::endl; 
   //std::cout << "diff\n" << ci_cjdag-cj_cidag.transpose() << std::endl; 
   //abort(); 

   //Mat Amat = Mat::Identity(n_site, n_site) - gf.transpose();  // A_{ij} = <c_i^dagger c_j>
   //std::cout << "nloc\n" << Amat.diagonal() << std::endl; 
   //
  /*local density*/
   //std::vector<double> nloc(n_site); 
   //for (unsigned i = 0; i< n_site; ++i){
   //    nloc[i] = (1.-gf(i,i)); 
   //}
   //measurements["nloc"] << nloc; 

   /*energies*/
   //double KinE = 0.0; 
   //BOOST_FOREACH(alps::graph_helper<>::site_descriptor const& sj, lattice.neighbors(si)) {//sj is neighbor of si 
   //   KinE += ci_cjdag(sj) + cj_cidag(sj); // assume hopping is one 
   //}
   //KinE *= n_cell; 
   //measurements["KinE"] << KinE;  

   double IntE = 0.; 
   BOOST_FOREACH(alps::graph_helper<>::site_descriptor const& sj, lattice.neighbors(si)) {
      IntE += std::pow(1.-ci_cjdag(si),2) - ci_cjdag(sj)*cj_cidag(sj);    //<n_i n_j>
   }
   IntE *= n_cell; 
   IntE = V*IntE - 0.25*V*n_bond; // (take care of constant shift )

   measurements["IntE"] << IntE; 
   //measurements["Energy"] << (KinE + IntE); 

   /*CDW structure factor*/   
   double delta; 
   double M2 = 0.; 
   for (site_t sj =0; sj< n_site; ++sj){  
        delta = (si == sj) ? 1.0 : 0.0;
        M2 += lattice.parity(si) * lattice.parity(sj) * 
              (std::pow(1. - ci_cjdag(si), 2) + (delta- cj_cidag(sj))*ci_cjdag(sj));
   }
   measurements["M2"] << 2.*M2/double(n_cell); 
}

