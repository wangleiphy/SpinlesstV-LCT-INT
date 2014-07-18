#include "interaction_expansion.hpp"
#include <ctime>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <boost/lexical_cast.hpp>
#include "buildK.h"

InteractionExpansion::InteractionExpansion(alps::params &parms, int node)
:alps::mcbase(parms,node),
Params(make_deprecated_parameters(parms)), 
lattice(Params),
max_order(boost::lexical_cast<unsigned int>(parms["MAX_ORDER"])),
n_site(lattice.num_sites()),
n_bond(lattice.num_bonds()),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned long)parms["THERMALIZATION"]),        
temperature(boost::lexical_cast<double>(parms["TEMPERATURE"])),                        
beta(1./temperature),  
V(boost::lexical_cast<double>(parms["V"])),                        
tlist(), 
vlist(), 
recalc_period(parms["RECALC_PERIOD"] | 500),
itime_max(boost::lexical_cast<itime_type>(parms["ITIME_MAX"])), 
nblock(boost::lexical_cast<itime_type>(parms["NBLOCKS"])),
steps_per_block(parms["STEPS_PER_BLOCK"] | 100),
blocksize(itime_max/nblock),
iblock(0),
direction(nblock==1? 0:1),
sweeps(0),
sign(1.), 
K_(buildK(lattice)),
gf(K_, beta, beta/boost::lexical_cast<double>(itime_max), nblock, blocksize, parms["UPDATE_REFRESH_PERIOD"] , parms["WRAP_REFRESH_PERIOD"])
//eng_(parms["SEED"] |42), 
//itime_rng(eng_, boost::uniform_int<itime_type>(0,itime_max)), 
//bond_rng(eng_, boost::uniform_int<site_type>(0,n_bond))
//random(eng_, boost::uniform_real<>()), 
{   
   
   //initialize ALPS observables
   initialize_observables();

   if(node==0) {
       print(std::cout); // print parameters to screen 
   }
}


void InteractionExpansion::update()// sweep in one block 
{

   for (itime_type i=0; i< nblock; ++i){

      sweeps++; // one sweep means try go through a block (with steps_per_block updates)
      interaction_expansion_step();                

      //if (tlist.size()== max_order)//TESTING: freeze the vertex configuration just do sweep 
      //    steps_per_block = 0; 

      iblock += direction; 

      //if hit the end, revert the sweep direction 
      if (iblock == nblock-1 || iblock == 0)
          direction *= -1; 

      //we jump to a new block and calculate gf at its time origin
      gf.wrap(iblock*blocksize, tlist, vlist); //this is necessary because otherwise we might jump over it_ some empty block 
 
      if(sweeps % recalc_period ==0)
         gf.rebuild(tlist, vlist);

    }
}


void InteractionExpansion::measure(){
  if (sweeps < therm_steps) 
   {
    //do nothing 
   }else{
    measure_observables();
   } 
}

double InteractionExpansion::fraction_completed() const {
    return (sweeps < therm_steps ? 0. : (sweeps - therm_steps)/double(nblock)/double(mc_steps));
}
