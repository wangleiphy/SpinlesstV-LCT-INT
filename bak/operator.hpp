#ifndef OPERATOR_HPP
#define OPERATOR_HPP

#include "types.h"

typedef class c_or_cdagger  //represents a creator operator at site s at time tau, site s has parity 
{ 
public:
  c_or_cdagger(const site_t s, const double parity, const itime_t t)
  :s_(s)
  ,parity_(parity)
  ,t_(t)
  {
  }

/*  
  const c_or_cdagger & operator=(const c_or_cdagger &c)
  {
    if(this != &c){
      z_ = c.z_;
      s_ = c.s_;
      t_ = c.t_;
    }
    return *this;
  }
  
  c_or_cdagger(const c_or_cdagger & c)
  {
    operator=(c);
  }
*/  

  inline const itime_t &t() const{return t_;}
  inline const double &parity() const {return parity_;}
  inline const site_t &s() const {return s_;}

//  inline const spin_t &flavor() const {return z_;}
//  inline spin_t &flavor() {return z_;}
//  inline itime_t &t() {return t_;}
//  inline void flavor(spin_t z){z_=z;}
//  inline void s(site_t s){s_=s;}
//  static void initialize_simulation(const alps::params &parms);
  
private:
  site_t s_;      //this operator site
  double parity_; // parity of this site 
  itime_t t_;     //its imaginary time point

} creator;

#endif
