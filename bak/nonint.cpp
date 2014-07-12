#include <Eigen/Dense>
#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/ngs.hpp>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/ngs/scheduler/parseargs.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "buildK.h"
#include "bgl.hpp"


double gf(const double E, const double beta, const double tau){ //<c(t) c^{+}(0)> with t>0 
         return (E>0.) ? exp(-E*tau)/(1.+exp(-beta*E)) :  exp((beta-tau)*E)/(1.+exp(beta*E)) ; // it is samething, to avoid overflow 
}


int main(int argc, char** argv){
   alps::parseargs options(argc, argv); 
   alps::params params(options.input_file);

   double beta = 1./boost::lexical_cast<double>(params["TEMPERATURE"]); 

   alps::Parameters Params(make_deprecated_parameters(params));
   alps::graph_helper<> lattice(Params); 
   unsigned int n_site(lattice.num_sites()); 

   //lattice helpers 
   typedef alps::graph_helper<>::site_descriptor site_descriptor; 
   typedef std::map<unsigned int, std::vector<site_descriptor> > DistanceMap;
   std::vector<DistanceMap> distmap(get_distmap(lattice)); 

   std::vector<unsigned int> shellsize(get_shellsize(distmap)); 

   Eigen::MatrixXi disttable(get_disttable(distmap, n_site)); 

   //hamiltonian 
   Eigen::MatrixXd K(buildK(lattice)); 

   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ces;
   ces.compute(K);
   
   Eigen::VectorXd v(n_site);
   for(unsigned int l=0; l<n_site; ++l) 
        v(l) = gf(ces.eigenvalues()(l), beta, 0.); 
   
  Eigen::MatrixXd g = ces.eigenvectors() * v.asDiagonal() * ces.eigenvectors().adjoint(); 

  unsigned int si =0; // the origin 
  std::vector<double> nncorr(shellsize.size());
  for (unsigned int sj= 0; sj< n_site; ++sj){
     unsigned int dist = disttable(si,sj); 
     double delta = (dist==0) ? 1.0 : 0.0;
     nncorr[dist] +=  (delta - g(sj,si)) *g(si, sj)/(double)shellsize[dist]; 
     std::cout << "dist, nncorr: " <<  dist << " " << (delta - g(sj,si)) *g(si, sj) << std::endl;  
  }

  for (unsigned int dist= 0; dist< nncorr.size(); ++dist){
      if (dist%2==1)
        std::cout << dist << " " << nncorr[dist] << std::endl ; 
  }

  return 0; 
}

