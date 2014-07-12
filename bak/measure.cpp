#include "interaction_expansion.hpp"
#include <vector>
#include <algorithm> 
//#include <boost/foreach.hpp>

void InteractionExpansion::measure_nncorrelation()
{

 unsigned int Msize = M.matrix().rows();

 double M2 = 0.; 
 double Kappa = 0.; 
 std::vector<double> nncorr(shellsize.size()); 

  for (unsigned int dist = 0; dist< shellsize.size(); ++dist){

           if (dist==0){//special treatment of si==sj, <(n-1/2)^2> = 0.25
                 nncorr[dist] += 0.25; 
                 M2 += 0.25;
                 Kappa += 0.25; 
              continue; 
           }

           site_t si = (unsigned int)(random()*n_site); 
           double tau = beta*random();

           DistanceMap dmap = distmap[si]; 
          
           unsigned int num_sites = dmap[dist].size();// number of sites in shell = dist 
           unsigned int j = random()*num_sites; 
           site_t sj = dmap[dist][j];  // pick up a random site in dist shell 
   
           std::vector<site_t> sites; 
           sites.push_back(si);
           sites.push_back(sj);
  
           Eigen::Matrix2d S = Eigen::Matrix2d::Zero(); // diagonal term is always zero ;  
            
           double parity = lattice.parity(sites[0]) * lattice.parity(sites[1]); 
           //compute the Green's functions with interpolation
           S(0,1) = bare_green_itime(0, sites[0], sites[1]); 
           //S(1,0) = bare_green_itime(0, sites[1], sites[0]); 
           S(1,0) = -S(0,1)* parity;   
           //for(unsigned int i=0; i<2; ++i){
           //    for(unsigned int j=0; j<2; ++j){
           //       S(i,j) = bare_green_itime(0, sites[i], sites[j]); 
           //    }
           //}
   
           if(Msize>0){
              Eigen::MatrixX2d Q(Msize,2);  // M*2
              Eigen::Matrix2Xd R(2,Msize);  // 2*M 
   
              for(unsigned int i=0; i<2; ++i){
               for(unsigned int j=0; j< Msize; ++j){
                    Q(j,i) = green0_spline(M.creators()[j].t()-tau, M.creators()[j].s(), sites[i]);
                    //R(i,j) = green0_spline(tau-M.creators()[j].t(), sites[i], M.creators()[j].s());
                    R(i,j) = -lattice.parity(sites[i]) * M.creators()[j].parity()* Q(j,i);//anti-symmetrization 
               }
             }
   
             S.noalias() -=  R * M.matrix() * Q; 
           }
           
           double ninj = S.determinant(); 
           //std::cout << "si, sj, tau, dist, ninj " << si << " " << sj  << " "<< tau << " " << dist << " " << ninj << std::endl; 
           //if (fabs(ninj)>10.){
           //   std::cout << "si, sj, tau, dist, ninj " << si << " " << sj  << " "<< tau << " " << dist << " " << ninj << std::endl; 
           //   std::cout << "creators: ";  
           //   for (unsigned int  i=0; i< M.creators().size(); ++i) {
           //     std::cout << M.creators()[i].s()<< "("<< M.creators()[i].t() << ")"  << ","; 
           //   }
           //   std::cout << std::endl; 
           //}

           nncorr[dist] += ninj;  //divid by number of sites in this shell 
           M2 += parity*ninj*num_sites;
           Kappa += ninj*num_sites;
  }
  
  //std::cout << "M2: "<<  2.*M2/(double)n_cell << std::endl; 

  measurements["nncorr"] << nncorr; 
  measurements["M2"] << 2.*M2/(double)n_cell; 
  measurements["Kappa"] << Kappa*beta; 
  measurements["IntE"] << V*nncorr[1]*(double)n_bond; // V(ni-1/2) (nj-1/2)
}

void InteractionExpansion::measure_local()
{

 unsigned int Msize = M.matrix().rows();

 std::vector<double> nloc(n_site); 
    
 for (unsigned int si=0; si< n_site; ++si){

          double tau = beta*random();
          double S = bare_green_itime(0, si, si); 

          if(Msize>0){
                Eigen::VectorXd Q(Msize);  
                Eigen::RowVectorXd R(Msize);
                
                 for(unsigned int j=0; j< Msize; ++j){
                      R(j) = green0_spline(tau-M.creators()[j].t(), si, M.creators()[j].s());
                      Q(j) = green0_spline(M.creators()[j].t()-tau, M.creators()[j].s(), si);
                }
                
                S -=  R * M.matrix() * Q; 
          }
          nloc[si] += S; 
 }

  std::ostream_iterator<double> out_it (std::cout,", ");
  std::cout << "nloc:\n";
  std::copy (nloc.begin(), nloc.end(), out_it);  
  std::cout << std::endl; 

}
