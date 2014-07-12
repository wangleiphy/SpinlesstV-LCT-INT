#ifndef BUILDK_H
#define BUILDK_H 

#include <boost/foreach.hpp>
#include <alps/lattice.h>
#include "types.h"

Eigen::MatrixXd buildK(const alps::graph_helper<>& lattice){


    //construct the hamiltonian then copy it to K 
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(lattice.num_sites(), lattice.num_sites()); 
    BOOST_FOREACH(const alps::graph_helper<>::site_descriptor& s1, lattice.sites()) {
        BOOST_FOREACH(alps::graph_helper<>::site_descriptor const& s2, lattice.neighbors(s1)) {
            K(s1, s2) = -1.0; 
        }
    }

   return K; 
}

#endif 
