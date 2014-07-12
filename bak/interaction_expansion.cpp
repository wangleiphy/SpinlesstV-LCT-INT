#include "interaction_expansion.hpp"
#include <ctime>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/model/hamiltonian_matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>
#include "buildK.h"
#include "bgl.hpp"

InteractionExpansion::InteractionExpansion(alps::params &parms, int node)
:alps::mcbase(parms,node),
Params(make_deprecated_parameters(parms)), 
lattice(Params),
max_order(boost::lexical_cast<unsigned int>(parms["MAX_ORDER"])),
//n_flavors(boost::lexical_cast<unsigned int>(parms["FLAVORS"])),
n_site(lattice.num_sites()),
n_bond(lattice.num_bonds()),
n_cell(n_site/2),
K_(buildK(lattice)), 
distmap(get_distmap(lattice)), 
disttable(get_disttable(distmap, n_site)), 
neighbors(get_neighbors(distmap, n_site, boost::lexical_cast<unsigned int>(parms["Nneighbors"]))),
shellsize(get_shellsize(distmap)), 
neighborshellsize(get_neighborshellsize(shellsize)), 
n_tau(boost::lexical_cast<unsigned int>(parms["N_TAU"])),
n_taumeasure(boost::lexical_cast<unsigned int>(parms["N_TAUMEASURE"])),
//n_self(parms["NSELF"] | boost::lexical_cast<unsigned int>(10*n_tau)),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned long)parms["THERMALIZATION"]),        
temperature(boost::lexical_cast<double>(parms["TEMPERATURE"])),                        
beta(1./temperature),  
timestepinv(n_tau*temperature),
timestep(1./timestepinv),
V(boost::lexical_cast<double>(parms["V"])),                        
//delta(boost::lexical_cast<double>(parms["DELTA"])),                        
recalc_period(parms["RECALC_PERIOD"] | 500),
measurement_period(parms["MEASUREMENT_PERIOD"] | 200),
M(), 
bare_green_itime(n_tau+1, n_site, K_, beta, timestep),
sweeps(0),
//weight(0.0),
sign(1.),
eta2(boost::lexical_cast<double>(parms["eta2"])),
eta4(boost::lexical_cast<double>(parms["eta4"])),
Zupdate(boost::lexical_cast<double>(parms["Zupdate"])),
ZtoW2(boost::lexical_cast<double>(parms["ZtoW2"])),
W2toZ(boost::lexical_cast<double>(parms["W2toZ"])),
W2toW4(boost::lexical_cast<double>(parms["W2toW4"])),
W4toW2(boost::lexical_cast<double>(parms["W4toW2"])),
ZtoW4(boost::lexical_cast<double>(parms["ZtoW4"])),
W4toZ(boost::lexical_cast<double>(parms["W4toZ"])),
Wupdate(boost::lexical_cast<double>(parms["Wupdate"])),
probs(),// empty vector 
coef2(2.*neighbors[0].size()/(double)n_cell), //constant ratio appear in worm create/destroy 
coef4(2.*std::pow(neighbors[0].size()/(double)n_cell, 3)),  
wormlength(0), // initialy we do not have worm
n_max(boost::lexical_cast<unsigned int>(parms["n_max"])), 
measure_unequaltime(boost::lexical_cast<bool>(parms["MEASURE_UNEQUALTIME"]) | false)
{
    probs.push_back(Zupdate);
    probs.push_back(2.*Zupdate);
    probs.push_back(2.*Zupdate+ZtoW2);
    probs.push_back(2.*Zupdate+ZtoW2+W2toZ);
    probs.push_back(2.*Zupdate+ZtoW2+W2toZ+W2toW4);
    probs.push_back(2.*Zupdate+ZtoW2+W2toZ+W2toW4+W4toW2);
    probs.push_back(2.*Zupdate+ZtoW2+W2toZ+W2toW4+W4toW2+ZtoW4);
    probs.push_back(2.*Zupdate+ZtoW2+W2toZ+W2toW4+W4toW2+ZtoW4+W4toZ);
    probs.push_back(2.*Zupdate+ZtoW2+W2toZ+W2toW4+W4toW2+ZtoW4+W4toZ+Wupdate);
    probs.push_back(2.*Zupdate+ZtoW2+W2toZ+W2toW4+W4toW2+ZtoW4+W4toZ+2.*Wupdate);

    //initialize ALPS observables
    initialize_observables();

   if(node==0) {
       print(std::cout); // print parameters to screen 
       update_params(parms); //write back a few generated params 
   }
}


void InteractionExpansion::update()
{
  for(unsigned int i=0;i<measurement_period;++i){
    sweeps++;
    interaction_expansion_step();                
    if(sweeps % recalc_period ==0)
      reset_perturbation_series();
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
