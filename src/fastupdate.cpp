#include "interaction_expansion.hpp"
#include <iterator> 

/*implement add vertex*/
double InteractionExpansion::add_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight)
{

  //Mat gtau = gf.G(itau, tlist, vlist) ; 

  Mat& gtau = gf.wrap(itau, tlist, vlist); // reference to its private member 
  std::cout << "gtau from wrap:\n"<< gtau<< std::endl; 
  std::cout << "gtau from scratch :\n"<< gf.G(itau, tlist, vlist)<< std::endl; 
  std::cout << "################################################ max diff:" <<  ((gtau - gf.G(itau, tlist, vlist)).cwiseAbs()).maxCoeff() << std::endl;

  site_type si = sites[0]; 
  site_type sj = sites[1]; 
  
  double gij = gf.U().row(si) * gtau * gf.Udag().col(sj); // rotate it to real space 

  double ratio = -4.* gij * gij; // gji = gij when they belongs to different sublattice 

  if(compute_only_weight){
       
  }else{

        //update gtau in eigen basis  
        gtau -= ( (gtau * gf.Udag().col(sj)) * (gf.U().row(si)* gtau) 
                 +(gtau * gf.Udag().col(si)) * (gf.U().row(sj)* gtau)
                )/gij;

        
        gtau -= 2.* (gtau* gf.Udag().col(si)) * gf.U().row(si) + 2.* (gtau*gf.Udag().col(sj)) * gf.U().row(sj); 

        //std::cout << "gtau from fastupdate:\n"<< gtau<< std::endl; 
   
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

 Mat& gtau = gf.wrap(itau, tlist, vlist); // reference to its private member, gtau is in eigen basis 
 std::cout << "gtau from wrap:\n"<< gtau<< std::endl; 
 std::cout << "gtau from scratch:\n"<<   gf.G(itau, tlist, vlist) << std::endl; 
 std::cout << "################################################ max diff:" <<  ((gtau - gf.G(itau, tlist, vlist)).cwiseAbs()).maxCoeff() << std::endl;
 
 site_type si = vlist[itau][0]; 
 site_type sj = vlist[itau][1]; 


 double gij = gf.U().row(si) * gtau * gf.Udag().col(sj); // rotate it to real space 

 double ratio = -4.* gij * gij; // gji = gij when they belongs to different sublattice 

 if(compute_only_weight){
    
 }else{

     //update gtau in eigen basis  
     gtau -= ( (gtau * gf.Udag().col(sj)) * (gf.U().row(si)* gtau) 
              +(gtau * gf.Udag().col(si)) * (gf.U().row(sj)* gtau)
             )/gij;

     // * U^\dagger V U  
     gtau -= 2.* (gtau* gf.Udag().col(si)) * gf.U().row(si) + 2.* (gtau*gf.Udag().col(sj)) * gf.U().row(sj); 

     //std::cout << "gtau from fastupdate:\n"<< gtau<< std::endl; 
     
     //update the vertex configuration 
     tlist.erase(itau); 
     vlist.erase(itau); 
     
     //std::cout << "gtau from scratch:\n"<< gf.G(itau, tlist, vlist) << std::endl; 
    
 }
     return ratio; 
}
