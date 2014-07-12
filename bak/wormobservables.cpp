#include "interaction_expansion.hpp"

///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables
void InteractionExpansion::initialize_observables() 
{
  measurements << alps::ngs::RealObservable("Sign")
               << alps::ngs::RealObservable("PertOrder")
               << alps::ngs::RealObservable("ZtoW2")
               << alps::ngs::RealObservable("W2toZ")
               << alps::ngs::RealObservable("W2toW4")
               << alps::ngs::RealObservable("W4toW2")
               << alps::ngs::RealObservable("ZtoW4")
               << alps::ngs::RealObservable("W4toZ")
               << alps::ngs::RealObservable("WormShift")
               ; 

 measurements  << alps::ngs::RealObservable("Z")
               << alps::ngs::RealObservable("W2")
               << alps::ngs::RealObservable("W4")
               << alps::ngs::RealObservable("IntE")
               << alps::ngs::RealObservable("Kappa")
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

   {
    std::stringstream obs_name;
    obs_name<<"WormAdd_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
   } 

   {
    std::stringstream obs_name;
    obs_name<<"WormRemoval_"<<i;
    measurements << alps::ngs::RealObservable(obs_name.str().c_str());
   }
 }


 for (unsigned int i=0; i< shellsize.size(); ++i){
   std::stringstream obs_name;
   obs_name<<"nncorr_"<<i;
   measurements << alps::ngs::RealObservable(obs_name.str().c_str());
 }

}

//this function is called whenever measurements should be performed.
void InteractionExpansion::measure_observables() 
{
  measurements["Sign"]<<sign;
  measurements["PertOrder"] << double(M.num_vertices()); // the true pert order 

  unsigned int dist = 0; 
  double nncorr = 0.; 

  if (wormlength ==0){
    measurements["Z"] << 1.;
    measurements["W2"] << 0.;
    measurements["W4"] << 0.;
    measurements["IntE"] << double(M.num_vertices());// the pert order measured in Z space 
    measurements["Kappa"] << 0.; 

  }else if (wormlength==2){

    measurements["Z"] << 0.;
    measurements["W2"] << 1.;
    measurements["W4"] << 0.;
    measurements["IntE"] << 0.; 

    //only here we have contribution to nncorr 
    site_t si = M.creators()[M.num_vertices()*2].s();
    site_t sj = M.creators()[M.num_vertices()*2+1].s(); 
    dist = disttable(si,sj); 
    nncorr = lattice.parity(si)*lattice.parity(sj)/(double)shellsize[dist]; 

    measurements["Kappa"] << lattice.parity(si)*lattice.parity(sj); 

  }else if (wormlength ==4){

    measurements["Z"] << 0.;
    measurements["W2"] << 0.;
    measurements["W4"] << 1.;
    measurements["IntE"] << 0.; 
    measurements["Kappa"] << 0.; 

  }

  for (unsigned int i=0; i< shellsize.size(); ++i){
    std::stringstream obs_name;
    obs_name<<"nncorr_"<<i;
    if (i==dist)
      measurements[obs_name.str().c_str()] << nncorr;
    else
      measurements[obs_name.str().c_str()] << 0.;
  }

}

//finial evaluation 
void InteractionExpansion::evaluate(results_type& results){

     results["IntE"] = (-1./beta)*results["IntE"]/results["Z"];
     results["Kappa"] = (0.25*beta*n_site)*results["Kappa"]/results["Z"]/eta2;

     results.insert("M2", results["W2"]/results["Z"]/eta2);
     results.insert("M4", results["W4"]/results["Z"]/eta4);

     results.insert("BinderRatio", (eta2*eta2/eta4)*results["W4"]*results["Z"]/(results["W2"]*results["W2"]));

     for (unsigned int i=0; i< shellsize.size(); ++i){
         std::stringstream obs_name;
         obs_name<<"nncorr_"<<i;
         results[obs_name.str().c_str()] *= 0.5*n_cell/eta2;  
         results[obs_name.str().c_str()] /= results["Z"];
     }
}
