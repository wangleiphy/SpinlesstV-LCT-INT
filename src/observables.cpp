#include "interaction_expansion.hpp"

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables
void InteractionExpansion::initialize_observables() 
{
  measurements << alps::ngs::RealObservable("Sign")
               << alps::ngs::RealObservable("PertOrder")
               << alps::ngs::RealObservable("Add")
               << alps::ngs::RealObservable("Removal")
               << alps::ngs::RealObservable("Shift")
               << alps::ngs::RealObservable("M2")
               << alps::ngs::RealObservable("IntE")
               << alps::ngs::RealObservable("KinE")
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
  measurements["Sign"]<<sign;
  measurements["PertOrder"] << double(tlist.size());

  measure_M2();

  //measure_vhist(); 

}

//finial evaluation 
void InteractionExpansion::evaluate(results_type& results){
    //empty 
}
