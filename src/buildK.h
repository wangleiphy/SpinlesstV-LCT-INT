#ifndef BUILDK_H
#define BUILDK_H 

#include <alps/lattice.h>
#include <boost/foreach.hpp>
#include "types.h"

///build the K matrix for the single particle hoppings
Mat buildK(const alps::graph_helper<>& lattice, const std::string BCmodifier){

    //construct the hamiltonian 
    Mat K = Mat::Zero(lattice.num_sites(), lattice.num_sites()); 
    
    alps::graph_helper<>::bond_iterator it, it_end;
    for (boost::tie(it, it_end) = lattice.bonds(); it != it_end; ++it) {
        alps::graph_helper<>::site_descriptor si = source(*it,lattice.graph()); 
        alps::graph_helper<>::site_descriptor sj = target(*it,lattice.graph()); 

        alps::graph_helper<>::boundary_crossing_type bc = get(alps::boundary_crossing_t(), lattice.graph(), *it);
        
        if (BCmodifier == "APBCX" )
          K(si,sj) = bc.crosses(0)==0 ? -1.0 : 1.0 ;//anti-periodic condition along x direction: if cross x, we revert sign of hopping  
        else
          K(si,sj) = -1.0; 

        K(sj,si) = K(si,sj);
    }

   //std::cout << "K:\n" << K << std::endl; 

   return K; 
}

Mat buildKtrial(const alps::graph_helper<>& lattice){

    //construct the hamiltonian 
    Mat K = Mat::Zero(lattice.num_sites(), lattice.num_sites()); 
    
    alps::graph_helper<>::bond_iterator it, it_end;
    for (boost::tie(it, it_end) = lattice.bonds(); it != it_end; ++it) {
        alps::graph_helper<>::site_descriptor si = source(*it,lattice.graph()); 
        alps::graph_helper<>::site_descriptor sj = target(*it,lattice.graph()); 

        alps::graph_helper<>::boundary_crossing_type bc = get(alps::boundary_crossing_t(), lattice.graph(), *it);
        
        //K(si,sj) = (bc.crosses(0)==0 && bc.crosses(1)==0) ? -1.0 : 0.0 ;//OBC condition along both directions   
        //K(si,sj) = (bc.crosses(0)==0 && bc.crosses(1)==0) ? -1.0 : 1.0 ;//APBC condition along both directions   
        K(si,sj) = (bc.crosses(0)==0) ? -1.0 : 1.0 ;//APBC condition along x directions   
        //K(si,sj) = -1.0; 
        K(sj,si) = K(si,sj);
    }

    /*
    //use my own random number generator because the one in class in not working right now 
    typedef boost::mt19937 engine_type;
    engine_type eng;
    boost::variate_generator<engine_type&, boost::uniform_real<> > random(eng, boost::uniform_real<>()); 

    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()){
         double delta  = 1E-3; 
         //double hopping = -1.0 + delta* (random() -0.5); // added random noise to break degeneracy
         double hopping = lattice.bond_type(b)==0 ? -1.-delta : -1.; 
         K(lattice.source(b), lattice.target(b)) = hopping; 
         K(lattice.target(b), lattice.source(b)) = hopping; 
    }

    */
   //std::cout << "K:\n" << K << std::endl; 

   return K; 
}

#endif 
