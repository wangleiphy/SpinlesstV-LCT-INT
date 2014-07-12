#include "interaction_expansion.hpp"
#include <iterator> 

/*implement add vertex*/
double InteractionExpansion::add_impl(const time_type tau, const std::vector<site_type>& sites, const bool compute_only_weight)
{

  Mat gtau = gf.G(tau, tlist, vlist); // gf at the current time 

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
  tlist.insert(tau); 
  vlist[tau] = sites; 

  return detratio; 
}


/*implement remove vertex in partition funciton sector*/
double InteractionExpansion::remove_impl(const unsigned vertex, const bool compute_only_weight)
{// this will only get call when there is a vertex 
 // vertices contains from small to large indices 

 tlist_type::const_iterator it(tlist.begin());    
 std::advance(it, vertex);
 double tau = *it; 

 Mat gtau = gf.G(tau, tlist, vlist); // gf at the current time 

 Eigen::Matrix2d S; 
 S(0, 0) = gtau(vlist[tau][0], vlist[tau][0]);
 S(0, 1) = gtau(vlist[tau][0], vlist[tau][1]);
 S(1, 0) = gtau(vlist[tau][1], vlist[tau][0]);
 S(1, 1) = gtau(vlist[tau][1], vlist[tau][1]);

 S = Eigen::Matrix2d::Identity() -2.0 *(Eigen::Matrix2d::Identity() - S); 
 double detratio = S.determinant(); 

 //return weight if we have nothing else to do
 if(compute_only_weight){
     return detratio; 
 }

 //update the vertex configuration 
 tlist.erase(tau); 
 vlist.erase(tau); 

 return detratio; 
 
}
