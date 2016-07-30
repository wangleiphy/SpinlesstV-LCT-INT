#include "interaction_expansion.hpp"

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables
void InteractionExpansion::initialize_observables() 
{
  measurements << alps::ngs::RealObservable("Sign")
               << alps::ngs::RealObservable("PertOrder")
               << alps::ngs::RealObservable("k")
               << alps::ngs::RealObservable("kPertOrder")
               << alps::ngs::RealObservable("kL")
               << alps::ngs::RealObservable("kR")
               << alps::ngs::RealObservable("kLkR")
               << alps::ngs::RealObservable("Add")
               << alps::ngs::RealObservable("Removal")
               << alps::ngs::RealObservable("M2")
               << alps::ngs::RealObservable("M2k")
               << alps::ngs::RealObservable("IntE")
               << alps::ngs::RealObservable("KinE")
               << alps::ngs::RealObservable("KinEk")
               << alps::ngs::RealObservable("Energy")
//               << alps::ngs::RealObservable("Kappa")
               ; 

  if (Add+Remove <1.0)
     measurements << alps::ngs::RealObservable("Shift"); 

  if (MEASURE_M4)
     measurements << alps::ngs::RealObservable("M4");

  //measurements << alps::ngs::RealVectorObservable("Vhist"); 

  measurements << alps::ngs::SimpleRealVectorObservable("nncorr"); 
  measurements << alps::ngs::RealObservable("Walltime"); 
}


//this function is called whenever measurements should be performed.
void InteractionExpansion::measure_observables() 
{

  double pert_order = tlist.size(); // number of vertices in the imaginary time 

  measurements["Sign"]<<sign;
  measurements["PertOrder"] << pert_order;
  
  //number of vertices in the center 
  tlist_type::const_iterator lower, upper; 
  lower = std::lower_bound (tlist.begin(), tlist.end(), window_lower); 
  upper = std::upper_bound (tlist.begin(), tlist.end(), window_upper, std::less_equal<itime_type>());  //equal is exclude

  double num_vertices = std::distance(lower, upper);  //number of vertices in the window

  measurements["k"] << num_vertices; 
  measurements["kPertOrder"] << num_vertices* pert_order; 

  //number of vertices on the left 
  lower = std::lower_bound (tlist.begin(), tlist.end(), 0); 
  upper = std::upper_bound (tlist.begin(), tlist.end(), itime_max/2, std::less_equal<itime_type>());  //equal is exclude

  double kL = std::distance(lower, upper); //number of vertices in the left part 
  double kR = pert_order - kL;  

  measurements["kL"] << kL;
  measurements["kR"] << kR;
  measurements["kLkR"] << kL*kR;

  measure_M2();
  if (MEASURE_M4){
      measure_M4(); 
  }

  //measure_vhist(); 
}

//finial evaluation 
void InteractionExpansion::evaluate(results_type& results){


    if (MEASURE_M4){
       results.insert("Binder", results["M4"]/(results["M2"]*results["M2"]));
    }

    results.insert("dlogKinEdV", (results["KinEk"]/results["KinE"] - results["PertOrder"])/V);
    results.insert("dlogM2dV", (results["M2k"]/results["M2"] - results["PertOrder"])/V);

    results.insert("IntE2", results["k"]/(-window_tau*n_site)); //Interaction energy per site 
    
    //second-order dervative of energy per site to V 
    results.insert("chiE", (results["kPertOrder"]- results["PertOrder"]*results["k"] - results["k"])/(-V*V*window_tau*n_site));

    results.insert("FS2", (results["kLkR"] - 0.25 * results["PertOrder"]*results["PertOrder"])/(V*V)); 
    results.insert("FS3", (results["kLkR"] - results["kL"]*results["kR"])/(V*V)); 

}
