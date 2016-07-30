#ifndef WEAK_COUPLING_H
#define WEAK_COUPLING_H

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/lattice.h>
#include <alps/alea.h>

#include <cmath>

#include "green_function.h"
//#include "green_function_Nmat.h"

class InteractionExpansion: public alps::mcbase
{
public:
  typedef boost::chrono::high_resolution_clock clock;

  InteractionExpansion(alps::params& p, int rank);
  ~InteractionExpansion() {}

  void update();
  void measure();
  double fraction_completed() const;

  //print progress 
  unsigned pertorder() const {return tlist.size();}
  unsigned long progress() const {return sweeps;}
  unsigned block() const {return iblock;}
  unsigned cycle() const {return cycles;}

  void evaluate(results_type& results);

  void test(); 
  void initialize_tvlist(); 

  using alps::mcbase::save;
  virtual void save(alps::hdf5::archive & ar) const;
  
  using alps::mcbase::load;
  virtual void load(alps::hdf5::archive & ar);

private:
  
  /*functions*/
  // in file io.cpp
  void print(std::ostream &os) const; //print parameters 
  
  /*the actual solver functions*/
  // in file solver.cpp
  void interaction_expansion_step();
  //void reset_perturbation_series();

  // in file model.cpp 
  // add and remove vertices 
  void add();
  void remove();
  void shift();

  // in file update.cpp:
  // add or remove vertex in partition funciton sector  
  double add_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight); 
  double remove_impl(const itime_type itau, const bool compute_only_weight);
  double shift_impl(const itime_type itau, const std::vector<site_type>& sites, const bool compute_only_weight);

  /*measurement functions*/
  // in file observables.cpp
  void initialize_observables();
  void measure_observables();
  // in file measure.cpp
  void measure_M2();
  void measure_M4();
  void measure_vhist(); 

  /*private member variables, constant throughout the simulation*/
  const alps::graph_helper<> lattice; 
  const unsigned int max_order;                        
    
  const unsigned n_site; 
  const unsigned n_bond; // number of *interaction* bond (fine when n.n. hopping and V )

  const boost::uint64_t mc_steps;                        
  const unsigned long therm_steps;                
  
  const time_type beta;  
  const double V;                        

  tlist_type tlist; //a list contains time where we have vertex 
  vlist_type vlist; //map from tau to sites 
  
  const unsigned recalc_period;                
  const unsigned measurement_period; 

  const itime_type itime_max;  
  const itime_type nblock; 
  const unsigned steps_per_block;        
  const itime_type blocksize;
  itime_type iblock; 
  int direction; 
  unsigned cycles; 
  unsigned long sweeps;        

  double sign;

  const time_type timestep; 
  const time_type window_tau;  // size of the window in which we perform measurement 

  const itime_type window_upper; //indices which indicate the start and end of the window 
  const itime_type window_lower; 


  Mat K_;    // the kinetic energy matrix 
  Mat Ktrial_; 
  Green_function gf; 

  //graph stuff 
  std::vector<DistanceMap> distmap;            //  vector<map from dist to vector<sites> >
  Eigen::MatrixXi          disttable;          //  table(si, sj) = dist  
  std::vector<unsigned>    shellsize;          //  number of sites in dist steps 


  const double Add; 
  const double Remove; 
  std::vector<double> probs; 
  const bool MEASURE_M4; 


  template<typename T>
  T randomint(const T i) {return random() * i;}//random int [0, i) 
    
};

#endif
