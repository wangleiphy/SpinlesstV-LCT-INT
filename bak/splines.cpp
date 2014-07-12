#include "interaction_expansion.hpp"

//as the name says: linearly interpolate between two points
//pass dxinv avoid division 
template<class X, class Y> inline Y linear_interpolate(const X x0, const X dxinv, const Y y0, const Y y1, const X x)
{
  //dxinv = 1/(x1-x0)
  return y0 + (x-x0)*dxinv*(y1-y0);
}

///Compute the Green's function G0 (the BARE) green's function between two points
double InteractionExpansion::green0_spline(const creator &cdagger, const creator &c) const
{
  //for the M matrix we need the bare temporal Green's function. 
  //For the dressed frequency GF we need the bare frequency GF.
  //here we get the bare temporal GF, interpolated, as needed for the M matrix.
  //we will receive this as an input into our solver later.
  itime_t delta_t = cdagger.t()-c.t();
  site_t site1 = cdagger.s();
  site_t site2 = c.s();

  //std::cout << "t1, t2, s1, s2 " << cdagger.t() << " " << c.t() << " " << site1 << " " << site2 << std::endl; 
  return green0_spline(delta_t, site1, site2);  
}


///Compute the bare green's function for a given site, and imaginary time.
double InteractionExpansion::green0_spline(const itime_t delta_t, const site_t site1, const site_t site2) const
{
  if(delta_t>=0.){
    int time_index = static_cast<int>(std::floor(delta_t*timestepinv));
    //std::cout << "detla_t, it1, it2 " << delta_t << " " << time_index_1 << " " << time_index_2 << std::endl; 
    return linear_interpolate(bare_green_itime.tau(time_index), timestepinv,
                              bare_green_itime(time_index,site1, site2),
                              bare_green_itime(time_index+1,site1, site2),delta_t);
  }else{
    //delta_t could be a small negative number, careful about precision issue 
    int time_index = static_cast<int>(std::floor(delta_t*timestepinv))+n_tau;
    //std::cout << "detla_t, it1, it2 " << delta_t << " " << time_index_1 << " " << time_index_2 << std::endl; 
    return -linear_interpolate(bare_green_itime.tau(time_index), timestepinv,
                               bare_green_itime(time_index,site1,site2),
                               bare_green_itime(time_index+1,site1,site2),delta_t+beta);
  }
}
