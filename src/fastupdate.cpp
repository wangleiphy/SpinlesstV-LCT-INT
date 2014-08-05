#include "interaction_expansion.hpp"
#include <iterator> 

/*implement add vertex*/
double InteractionExpansion::add_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight)
{
    
  gf.wrap(itau, tlist, vlist);  

  site_type si = sites[0]; 
  site_type sj = sites[1]; 
  
  double gij = gf.gij(si, sj); // rotate it to real space 

  double ratio = -4.* gij * gij; // gji = gij when they belongs to different sublattice 

  if(compute_only_weight){
       
  }else{

        gf.update(si, sj, gij, gij); 
        //std::cout << "gtau from fastupdate:\n"<< gf.gtau() << std::endl; 
   
        //update the vertex configuration 
        tlist.insert(itau); 
        vlist[itau] = sites; 
   
        //std::cout << "gtau from scratch:\n"<<  gf.G(itau, tlist, vlist) << std::endl; 
  }

  return ratio; 
}


/*implement remove vertex at time itau*/
double InteractionExpansion::remove_impl(const itime_type itau, const bool compute_only_weight)
{// this will only get call when there is a vertex 

 gf.wrap(itau, tlist, vlist);  

 site_type si = vlist[itau][0]; 
 site_type sj = vlist[itau][1]; 

 double gij = gf.gij(si, sj);  

 double ratio = -4.* gij * gij; // gji = gij when they belongs to different sublattice 

 if(compute_only_weight){
    
 }else{
        
     if (tlist.size()> 1){
         gf.update(si, sj, gij, gij); 
     }else{ //since we will get a empty list as use this opputunity to reset all memory 
         //std::cout << "empty vertex occur" << std::endl; 
         gf.init_without_vertex(); 
     }

     //std::cout << "gtau from fastupdate:\n"<< gf.gtau() << std::endl; 
     
     //update the vertex configuration 
     tlist.erase(itau); 
     vlist.erase(itau); 
     //std::cout << "gtau from scratch:\n"<< gf.G(itau, tlist, vlist) << std::endl; 
    
 }
     return ratio; 
}

/*implement shift site index in a vertex*/
double InteractionExpansion::shift_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight)
{

  gf.wrap(itau, tlist, vlist);  
    
  site_type si = sites[0]; 
  site_type sj = sites[1]; 
  site_type sjprime = sites[2]; 
  
  double gjjprime = gf.gij(sj, sjprime); // rotate it to real space 

  double ratio = 4.* gjjprime * gjjprime; // gji = -gij when they belongs to same sublattice 

  if(compute_only_weight){
       
  }else{

        gf.update(sj, sjprime, gjjprime, -gjjprime); 
        //std::cout << "gtau from fastupdate:\n"<< gf.gtau() << std::endl; 
   
        //update the vertex configuration 
        std::vector<site_type> new_sites; 
        new_sites.push_back(si); 
        new_sites.push_back(sjprime); 

        vlist[gf.itau()] = new_sites; 
        //std::cout << "gtau from scratch:\n"<<  gf.G(itau, tlist, vlist) << std::endl; 
  }

  return ratio; 
}
