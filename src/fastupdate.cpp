#include "interaction_expansion.hpp"
#include <iterator> 

/*implement add vertex*/
double InteractionExpansion::add_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight)
{
    

  Mat gtau = gf.G(itau, tlist, vlist); // gf at the current time 

  //std::cout << "gtau:\n"<< gtau<< std::endl; 

  Eigen::Matrix2d S; 
  S(0, 0) = gtau(sites[0], sites[0]);
  S(0, 1) = gtau(sites[0], sites[1]);
  S(1, 0) = gtau(sites[1], sites[0]); 
  S(1, 1) = gtau(sites[1], sites[1]);

  S = Eigen::Matrix2d::Identity() -2.0 *(Eigen::Matrix2d::Identity() - S); 

  double detratio = S.determinant(); 

  //return weight if we have nothing else to do
  if(compute_only_weight){
      return detratio; 
  }

  //update the vertex configuration 
  tlist.insert(itau); 
  vlist[itau] = sites; 

  return detratio; 
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

 Eigen::Matrix2d S; 
 S(0, 0) = gtau(si, si);
 S(0, 1) = gtau(si, sj);
 S(1, 0) = gtau(sj, si);
 S(1, 1) = gtau(sj, sj);

 S = Eigen::Matrix2d::Identity() -2.0 *(Eigen::Matrix2d::Identity() - S); 
 double detratio = S.determinant(); 

 //return weight if we have nothing else to do
 if(compute_only_weight){
     return detratio; 
 }

 //update the vertex configuration 
 tlist.erase(itau); 
 vlist.erase(itau); 

 return detratio; 
 
}
