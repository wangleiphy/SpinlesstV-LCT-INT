#include <alps/ngs.hpp>
#include <alps/ngs/scheduler/parseargs.hpp>
#include "interaction_expansion.hpp"

void InteractionExpansion::test(){
    
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
                std::cout << "add vertex with detratio: " << detratio<< std::endl; 
                std::cout << "number of vertices: " << tlist.size() << std::endl; 
     
}


int main(int argc, char** argv){
   alps::parseargs options(argc, argv); 
   alps::params params(options.input_file);

   unsigned Ntry = 1; 
   for (unsigned itry = 0; itry < Ntry ; ++itry) {
        InteractionExpansion sim(params, 0); 
        sim.test(); 
   }

   return 0; 
}
