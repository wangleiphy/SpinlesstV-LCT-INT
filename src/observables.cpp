#include "interaction_expansion.hpp"

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables
void InteractionExpansion::initialize_observables() 
{
  measurements << alps::ngs::RealObservable("Sign")
               << alps::ngs::RealObservable("PertOrder")
               << alps::ngs::RealObservable("k")
               << alps::ngs::RealObservable("kPertOrder")
               << alps::ngs::RealObservable("Add")
               << alps::ngs::RealObservable("Removal")
               << alps::ngs::RealObservable("Shift")
               << alps::ngs::RealObservable("M2")
               << alps::ngs::RealObservable("M2k")
               << alps::ngs::RealObservable("IntE")
               << alps::ngs::RealObservable("KinE")
               << alps::ngs::RealObservable("KinEk")
               << alps::ngs::RealObservable("Energy")
//               << alps::ngs::RealObservable("Kappa")
               ; 

  //measurements << alps::ngs::RealVectorObservable("Vhist"); 

  measurements << alps::ngs::RealVectorObservable("nncorr"); 
  measurements << alps::ngs::RealObservable("Walltime"); 
}


//this function is called whenever measurements should be performed.
void InteractionExpansion::measure_observables() 
{

  double pert_order = tlist.size(); // number of vertices in the imaginary time 

  measurements["Sign"]<<sign;
  measurements["PertOrder"] << pert_order;

  tlist_type::const_iterator lower, upper; 
  lower = std::lower_bound (tlist.begin(), tlist.end(), window_lower); 
  upper = std::upper_bound (tlist.begin(), tlist.end(), window_upper, std::less_equal<itime_type>());  //equal is exclude

  double num_vertices = std::distance(lower, upper);  //number of vertices in the window

  measurements["k"] << num_vertices; 
  measurements["kPertOrder"] << num_vertices* pert_order; 

  measure_M2();

  //measure_vhist(); 

}

//finial evaluation 
void InteractionExpansion::evaluate(results_type& results){
    results.insert("dlogKinEdV", (results["KinEk"]/results["KinE"] - results["PertOrder"])/V);
    results.insert("dlogM2dV", (results["M2k"]/results["M2"] - results["PertOrder"])/V);

    results.insert("dEnergydV2", (results["kPertOrder"]- results["PertOrder"]*results["k"] - results["k"])/(-V*V*n_bond *window_tau));

}
