#include "interaction_expansion.hpp"

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables
void InteractionExpansion::initialize_observables() 
{
  measurements << alps::ngs::RealObservable("Sign")
               << alps::ngs::RealObservable("PertOrder")
               ; 

 for (unsigned int i=1; i<= n_max; ++i){
   {
    std::stringstream obs_name;
    obs_name<<"VertexAdd_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
   }

   {
    std::stringstream obs_name;
    obs_name<<"VertexRemoval_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
   }
 }

    measurements << alps::ngs::RealObservable("M2")
                 << alps::ngs::RealObservable("Kappa")
                 << alps::ngs::RealObservable("IntE");

    measurements << alps::ngs::RealVectorObservable("nncorr"); 

    if (measure_unequaltime)
        measurements << alps::ngs::RealVectorObservable("gf");
    //               << alps::ngs::RealVectorObservable("ntaun");
}


//this function is called whenever measurements should be performed.
void InteractionExpansion::measure_observables() 
{
  measurements["Sign"]<<sign;
  measurements["PertOrder"] << double(M.num_vertices());

  //measure_local();  local density is always 0.5 because RMQ is anti-symmetric 
  measure_nncorrelation(); 
  //measure_equaltime(); 

  if (measure_unequaltime)
    measure_gf(); 
  //measure_ntaun(); 

}

//finial evaluation 
void InteractionExpansion::evaluate(results_type& results){
    //empty 
}
