#include <alps/ngs.hpp>
#include <alps/ngs/scheduler/parseargs.hpp>
#include "interaction_expansion.hpp"

void InteractionExpansion::test(){
    
      unsigned Ntry = 2; 

      for (unsigned itry = 0; itry < Ntry ; ++itry) {

            //if (random() < 0.5 ) {//adds ONE vertex
            //if (itry%2==0 ) {//adds ONE vertex
            if ( true ) {//adds ONE vertex
                std::vector<site_type> sites;  
            
                alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
                sites.push_back(lattice.source(b));
                sites.push_back(lattice.target(b));
            
                itime_type itau = randomint(itime_max);
            
                std::cout << "#######################"  << std::endl; 
                std::cout << "weight before: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 

                double detratio = add_impl(itau, sites, false);  

                std::cout << "weight after: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 
                std::cout << "itau, b:"  << itau << " " << b << std::endl; 
                std::cout << "tlist: "; 
                std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
                std::cout << std::endl; 
                std::cout << "itry, add vertex with detratio: " << itry << " " << detratio<< std::endl; 
                std::cout << "number of vertices: " << tlist.size() << std::endl; 
            
            }else{
                if (tlist.size()< 1) continue; 
                unsigned vertex = randomint(tlist.size());

                tlist_type::const_iterator it = tlist.begin();  
                std::advance(it, vertex);  
                itime_type itau = *it; 

                std::cout << "#######################"  << std::endl; 
                std::cout << "weight before: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 

                double detratio = remove_impl(itau, false); 
                
                std::cout << "weight after: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 
                std::cout << "v:" << vertex  << std::endl; 
                std::cout << "tlist: "; 
                std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
                std::cout << std::endl; 
                std::cout << "itry, remove vertex with detratio: " << itry << " "<<  detratio<< std::endl; 
                std::cout << "number of vertices: " << tlist.size() << std::endl; 
            }
      }
}


int main(int argc, char** argv){
   alps::parseargs options(argc, argv); 
   alps::params params(options.input_file);

   InteractionExpansion sim(params, 0); 
   std::cout << "initialization done" << std::endl; 
   sim.test(); 
   return 0; 
}

