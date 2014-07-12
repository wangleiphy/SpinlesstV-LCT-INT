#ifndef WEAK_COUPLING_H
#define WEAK_COUPLING_H

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/lattice.h>
#include <alps/alea.h>

#include <cmath>
#include "green_function.h"
#include "types.h"
#include "mmatrix.hpp"
#include "operator.hpp"
#include <boost/chrono.hpp>

/*types*/
class c_or_cdagger;

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
  unsigned int pertorder() const {return M.num_vertices();}; 
  unsigned long progress() const {return sweeps;};        

  void build_matrix(); 
  void test(); 

  void evaluate(results_type& results);


private:
  
  /*functions*/
  // in file io.cpp
  void print(std::ostream &os) const; //print parameters 
  void update_params(alps::params &parms) const; 
  
  /*green's function*/
  // in file spines.cpp
  double green0_spline(const creator &cdagger, const creator &c) const;
  double green0_spline(const itime_t delta_t, const site_t site1, const site_t site2) const;
  
  /*the actual solver functions*/
  // in file solver.cpp
  void interaction_expansion_step();
  void reset_perturbation_series();

  // in file model.cpp 
  // add and remove vertices 
  void add();
  void remove();

  // in file worm.cpp 
  // create and destroy worm, it forms a triangle between Z, M2 and M4 
  void Z_to_W2(); 
  void W2_to_Z(); 

  void W2_to_W4(); 
  void W4_to_W2(); 

  void Z_to_W4(); 
  void W4_to_Z(); 

  //add and remove vertex in worm space 
  void wormadd(); 
  void wormremove(); 
  //shift worm 
  void wormshift(); 
  //create/destroy worm by open a vertex 
  //void wormopen(); 
  //void wormclose(); 
 

  // in file update.cpp:
  // add or remove vertex in partition funciton sector  
  double add_impl(const std::vector<double>& taus, const std::vector<site_t>& sites, const bool compute_only_weight); 
  double remove_impl(const std::vector<unsigned int>& vertices, const bool compute_only_weight);

  //in file wormupdate.cpp:
  //create or destroy the worm
  double wormcreate_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight);
  double wormdestroy_impl(const bool compute_only_weight, const unsigned int Morder);
  
  //add or remove vertex in the presence of worm
  double wormadd_impl(const std::vector<double>& taus, const std::vector<site_t>& sites, const bool compute_only_weight);
  double wormremove_impl(const std::vector<unsigned int>& vertices, const bool compute_only_weight);
  //shift a worm 
  double wormshift_impl(const double tau, const std::vector<site_t>& sites, const bool compute_only_weight);
  //open/close a vertex to get worm  
  double wormopen_impl(const unsigned int vertex_nr, const bool compute_only_weight); 
  double wormclose_impl(const bool compute_only_weight);
  
  /*measurement functions*/
  // in file observables.cpp
  void initialize_observables();
  void measure_observables();
  //in file measure.cpp
  void measure_local(); 
  void measure_nncorrelation(); 
  //in file unequaltime.cpp 
  void measure_gf();     // <c(t) c^{+}(0)>
  void measure_ntaun();  // <(n(t)-1/2)(n(0)-1/2)>
 
  /*private member variables, constant throughout the simulation*/
  const alps::Parameters Params;
  const alps::graph_helper<> lattice; 
  const unsigned int max_order;                        
  //const spin_t n_flavors;                          //number of flavors 
  const site_t n_site;                               //number of sites
  const site_t n_bond;                               //number of *interaction* bond (fine when n.n. hopping and V )
  const site_t n_cell;                               //number of unit cells = n_site/2 
  Eigen::MatrixXd K_;  // the kinetic energy matrix 
  
  //graph stuff 
  std::vector<DistanceMap> distmap;             //  vector<map(dist:vector<sites>)>
  Eigen::MatrixXi          disttable;           //  table(si, sj) = dist  
  std::vector<std::vector<site_t> > neighbors;  //neighbor list for each site within Nneighbors hoppings 
  std::vector<unsigned int> shellsize;          //number of sites in dist steps 
  std::vector<unsigned int> neighborshellsize;  //number of sites in dist+1 and dist-1 steps

//  const frequency_t n_matsubara;        //number of matsubara freq
//  const frequency_t n_matsubara_measurements;        //number of measured matsubara freq
  const itime_index_t n_tau;                        //number of imag time slices
  const itime_index_t n_taumeasure;                 // number of imag time where we do measurement 
//  const frequency_t n_self;                        //number of self energy (W) binning points
  const boost::uint64_t mc_steps;                        
  const unsigned long therm_steps;                
  
  const itime_t temperature;                               
  const itime_t beta;  
  const itime_t timestepinv;// n_tau * temperature 
  const itime_t timestep;   // 1/(n_tau * temperature) 
  const double V;                        
  //const double delta; 
  
  const unsigned int recalc_period;                
  const unsigned int measurement_period;        
  //const unsigned int convergence_check_period;        
  
  //important M matrix  
  m_matrix M;

  /*private member variables*/
  itime_green_function_t bare_green_itime;
    
  unsigned long sweeps;        
  //double weight;
  double sign;

  const double eta2; //coef before M2
  const double eta4; //coef before M4 
  const double Zupdate; 
  const double ZtoW2;
  const double W2toZ; 
  const double W2toW4;
  const double W4toW2;
  const double ZtoW4;
  const double W4toZ;
  const double Wupdate; //(Zupdate + Wupdate ) *2 + M2create + M2destroy + M4create + M4destroy + Wshift = 1
  std::vector<double> probs; 

  const double coef2; 
  const double coef4; 

  unsigned int wormlength; 

  const unsigned int n_max; 
  const bool measure_unequaltime; 

  unsigned int randomint(const unsigned int i) {return random() * i;}//random int [0, i) 

};

#endif
