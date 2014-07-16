#include <alps/ngs.hpp>
#include <alps/ngs/scheduler/parseargs.hpp>
#include "interaction_expansion.hpp"

void InteractionExpansion::test(){
    
      unsigned Ntry_perblock = 2; 
        
      for (unsigned b = 0; b < nblock; ++b ){// sweep through blocks 
        for (unsigned itry = 0; itry < Ntry_perblock ; ++itry) {

            //if (random() < 0.5 ) {//adds ONE vertex
            if (itry%2==0 ) {//adds ONE vertex
            //if ( true ) {//adds ONE vertex
                std::vector<site_type> sites;  
            
                alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
                sites.push_back(lattice.source(b));
                sites.push_back(lattice.target(b));
            
                itime_type itau = iblock*blocksize + randomint(blocksize);// a random time inside this block 
            
                std::cout << "######add#################"  << std::endl; 
                std::cout << "weight before: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 
                std::cout << "itau, b:"  << itau << " " << b << std::endl; 

                double detratio = add_impl(itau, sites, false);  

                std::cout << "weight after: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 
                std::cout << "tlist: "; 
                std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
                std::cout << std::endl; 
                std::cout << "itry, add vertex with detratio: " << itry << " " << std::setprecision(9)  << detratio<< std::endl; 
                std::cout << "number of vertices: " << tlist.size() << std::endl; 
            
            }else{
                if (tlist.size()< 1) continue; 

                tlist_type::const_iterator lower, upper; 
                lower = std::lower_bound (tlist.begin(), tlist.end(), iblock*blocksize); 
                upper = std::upper_bound (tlist.begin(), tlist.end(), (iblock+1)*blocksize, std::less_equal<itime_type>());  //equal is exclude
   
                unsigned num_vertices = std::distance(lower, upper); //number of vertices in this block
   
                if(num_vertices < 1){
                    return; 
                }
                
                std::advance(lower, randomint(num_vertices)); //the vertex to remove 
                itime_type itau = *lower; 

                std::cout << "####remove###################"  << std::endl; 
                std::cout << "weight before: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 

                double detratio = remove_impl(itau, false); 
                
                std::cout << "weight after: " << 1./gf.G(0, tlist, vlist).determinant() << std::endl; 
                std::cout << "tlist: "; 
                std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
                std::cout << std::endl; 
                std::cout << "itry, remove vertex with detratio: " << itry << " "<<  std::setprecision(9)  << detratio<< std::endl; 
                std::cout << "number of vertices: " << tlist.size() << std::endl; 
            }
        }   

            gf.rebuild(tlist, vlist);
            
            iblock += direction; 
            //we jump to a new block and calculate gf at its time origin
            gf.wrap(iblock*blocksize, tlist, vlist); //this is necessary because otherwise we might jump over it_ some empty block 
            //if hit the end, revert the sweep direction 
            if (iblock == nblock-1 || iblock == 0)
                direction *= -1; 
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

