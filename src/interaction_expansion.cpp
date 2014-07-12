#include "interaction_expansion.hpp"
#include <ctime>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>
#include "buildK.h"

InteractionExpansion::InteractionExpansion(alps::params &parms, int node)
:alps::mcbase(parms,node),
Params(make_deprecated_parameters(parms)), 
lattice(Params),
max_order(boost::lexical_cast<unsigned int>(parms["MAX_ORDER"])),
n_bond(lattice.num_bonds()),
K_(buildK(lattice)), 
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned long)parms["THERMALIZATION"]),        
temperature(boost::lexical_cast<double>(parms["TEMPERATURE"])),                        
beta(1./temperature),  
V(boost::lexical_cast<double>(parms["V"])),                        
tlist(), 
vlist(), 
gf(K_, beta), 
recalc_period(parms["RECALC_PERIOD"] | 500),
measurement_period(parms["MEASUREMENT_PERIOD"] | 200),
sweeps(0),
sign(1.)
{
    //initialize ALPS observables
    initialize_observables();

   if(node==0) {
       print(std::cout); // print parameters to screen 
   }
}


void InteractionExpansion::update()
{
  for(unsigned int i=0;i<measurement_period;++i){
    sweeps++;
    interaction_expansion_step();                
    //if(sweeps % recalc_period ==0)
    //  reset_perturbation_series();
  }
}

void InteractionExpansion::measure(){
  if (sweeps  <  therm_steps) 
   {
    //do nothing 
   }else{
    measure_observables();
   } 
}

double InteractionExpansion::fraction_completed() const {
    return (sweeps < therm_steps ? 0. : ( sweeps - therm_steps )/double(measurement_period)/ double(mc_steps));
}
