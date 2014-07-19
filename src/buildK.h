#ifndef BUILDK_H
#define BUILDK_H 

#include <alps/lattice.h>
#include <boost/foreach.hpp>
#include "types.h"

Mat buildK(const alps::graph_helper<>& lattice, const std::string BC){

    //use my own random number generator because the one in class in not working right now 
    typedef boost::mt19937 engine_type;
    engine_type eng;
    boost::variate_generator<engine_type&, boost::uniform_real<> > random(eng, boost::uniform_real<>()); 

    //construct the hamiltonian 
    Mat K = Mat::Zero(lattice.num_sites(), lattice.num_sites()); 
    
    alps::graph_helper<>::bond_iterator it, it_end;
    for (boost::tie(it, it_end) = lattice.bonds(); it != it_end; ++it) {
        alps::graph_helper<>::site_descriptor si = source(*it,lattice.graph()); 
        alps::graph_helper<>::site_descriptor sj = target(*it,lattice.graph()); 

        alps::graph_helper<>::boundary_crossing_type bc = get(alps::boundary_crossing_t(), lattice.graph(), *it);
        
        if (BC == "APBCX" )
          K(si,sj) = bc.crosses(0)==0 ? -1.0 : 1.0 ;//anti-periodic condition along x direction: if cross x, we revert sign of hopping  
        else
          K(si,sj) = -1.0; 

        K(sj,si) = K(si,sj);
    }

   /*
    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()){
         double hopping = -1.0; // + 0.0001 * (random() -0.5); // added random noise to break degeneracy
         K(lattice.source(b), lattice.target(b)) = hopping; 
         K(lattice.target(b), lattice.source(b)) = hopping; 
    }
   */


   //std::cout << "K:\n" << K << std::endl; 

   return K; 
}

#endif 
