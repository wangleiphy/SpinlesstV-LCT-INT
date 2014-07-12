#ifndef WEAK_COUPLING_H
#define WEAK_COUPLING_H

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/lattice.h>
#include <alps/alea.h>

#include <cmath>
#include "green_function.h"

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
  unsigned pertorder() const {return tlist.size();}; 
  unsigned long progress() const {return sweeps;};        
  void evaluate(results_type& results);


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

  // in file update.cpp:
  // add or remove vertex in partition funciton sector  
  double add_impl(const double tau, const std::vector<site_type>& sites, const bool compute_only_weight); 
  double remove_impl(const unsigned vertex, const bool compute_only_weight);

  /*measurement functions*/
  // in file observables.cpp
  void initialize_observables();
  void measure_observables();

  /*private member variables, constant throughout the simulation*/
  const alps::Parameters Params;
  const alps::graph_helper<> lattice; 
  const unsigned int max_order;                        

  const unsigned n_bond; // number of *interaction* bond (fine when n.n. hopping and V )
  Eigen::MatrixXd K_;    // the kinetic energy matrix 

  const boost::uint64_t mc_steps;                        
  const unsigned long therm_steps;                
  
  const time_type temperature;                               
  const time_type beta;  
  const double V;                        
  tlist_type tlist; //a list contains time where we have vertex 
  vlist_type vlist; //map from tau to sites 
  Green_function gf; 
  
  const unsigned recalc_period;                
  const unsigned measurement_period;        
  unsigned long sweeps;        

  double sign;

  unsigned int randomint(const unsigned int i) {return random() * i;}//random int [0, i) 

};

#endif
