#include "interaction_expansion.hpp"
#include <iterator> 

/*implement add vertex*/
double InteractionExpansion::add_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight)
{
    
  Mat gtau = gf.G(itau, tlist, vlist); // gf at the current time 
  //std::cout << "gtau:\n"<< gtau<< std::endl; 

  site_type si = sites[0]; 
  site_type sj = sites[1]; 
  double ratio = -4.* gtau(si, sj) *  gtau(sj, si); 

  if(compute_only_weight){
       
  }else{

        //update gtau 
        gtau -= gtau.col(sj)*gtau.row(si)/gtau(si, sj) + gtau.col(si) * gtau.row(sj)/gtau(sj, si);

        gtau.col(si) *= -1.; 
        gtau.col(sj) *= -1.; 

        std::cout << "gtau from fastupdate:\n"<< gtau<< std::endl; 
   
        //update the vertex configuration 
        tlist.insert(itau); 
        vlist[itau] = sites; 
   
        gtau = gf.G(itau, tlist, vlist);
        std::cout << "gtau from scratch:\n"<< gtau<< std::endl; 
        
   
  }

  return ratio; 
}


/*implement remove vertex in partition funciton sector*/
double InteractionExpansion::remove_impl(const unsigned vertex, const bool compute_only_weight)
{// this will only get call when there is a vertex 
 // vertices contains from small to large indices 

 tlist_type::const_iterator it(tlist.begin());    
 std::advance(it, vertex);
 itime_type itau = *it; 

 Mat gtau = gf.G(itau, tlist, vlist); // gf at the current time 
 //std::cout << "gtau:\n"<< gtau<< std::endl; 
 
 site_type si = vlist[itau][0]; 
 site_type sj = vlist[itau][1]; 

 double ratio = -4.* gtau(si, sj) *  gtau(sj, si); 

 if(compute_only_weight){
    
 }else{

     //update gtau 
     gtau -= gtau.col(sj)*gtau.row(si)/gtau(si, sj) + gtau.col(si) * gtau.row(sj)/gtau(sj, si);

     gtau.col(si) *= -1.; 
     gtau.col(sj) *= -1.; 
     std::cout << "gtau from fastupdate:\n"<< gtau<< std::endl; 
     
     //update the vertex configuration 
     tlist.erase(itau); 
     vlist.erase(itau); 
     
     gtau = gf.G(itau, tlist, vlist);
     std::cout << "gtau from scratch:\n"<< gtau<< std::endl; 
    
 }
     return ratio; 
}
