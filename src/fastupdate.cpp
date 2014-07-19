#include "interaction_expansion.hpp"
#include <iterator> 

/*implement add vertex*/
double InteractionExpansion::add_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight)
{

  gf.wrap(itau, tlist, vlist);  
    
  //const Mat gtau = gf.gtau(); 
  //std::cout << "gtau from wrap:\n" << gtau << std::endl; 
  //std::cout << "gtau from scratch :\n"<< gf.G(itau, tlist, vlist)<< std::endl; 
  //std::cout << "################################################ max diff:" <<  ((gtau- gf.G(itau, tlist, vlist)).cwiseAbs()).maxCoeff() << std::endl;

  site_type si = sites[0]; 
  site_type sj = sites[1]; 
  
  double gij = gf.gij(si, sj); // rotate it to real space 

  double ratio = -4.* gij * gij; // gji = gij when they belongs to different sublattice 

  if(compute_only_weight){
       
  }else{

        gf.update(si, sj, gij); 
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
 // vertices contains from small to large indices 


 //Mat gtau = gf.G(itau, tlist, vlist) ; 
 gf.wrap(itau, tlist, vlist); // reference to its private member, gtau is in eigen basis 


 //const Mat gtau = gf.gtau(); 
 //std::cout << "gtau from wrap:\n" << gtau << std::endl; 
 //std::cout << "gtau from scratch :\n"<< gf.G(itau, tlist, vlist)<< std::endl; 
 //std::cout << "################################################ max diff:" <<  (( gtau - gf.G(itau, tlist, vlist)).cwiseAbs()).maxCoeff() << std::endl;

 site_type si = vlist[itau][0]; 
 site_type sj = vlist[itau][1]; 

 double gij = gf.gij(si, sj);  

 double ratio = -4.* gij * gij; // gji = gij when they belongs to different sublattice 

 if(compute_only_weight){
    
 }else{
        
     if (tlist.size()> 1)
         gf.update(si, sj, gij); 
     else //since we will get a empty list as use this opputunity to reset all memory 
         gf.init_without_vertex(); 

     //std::cout << "gtau from fastupdate:\n"<< gf.gtau() << std::endl; 
     
     //update the vertex configuration 
     tlist.erase(itau); 
     vlist.erase(itau); 
     //std::cout << "gtau from scratch:\n"<< gf.G(itau, tlist, vlist) << std::endl; 
    
 }
     return ratio; 
}
