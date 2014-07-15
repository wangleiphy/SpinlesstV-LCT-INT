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
n_site(lattice.num_sites()),
n_bond(lattice.num_bonds()),
K_(buildK(lattice)), 
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned long)parms["THERMALIZATION"]),        
temperature(boost::lexical_cast<double>(parms["TEMPERATURE"])),                        
beta(1./temperature),  
V(boost::lexical_cast<double>(parms["V"])),                        
tlist(), 
vlist(), 
recalc_period(parms["RECALC_PERIOD"] | 500),
nblock(parms["NBLOCKS"] | 1),
steps_per_block(parms["STEPS_PER_BLOCK"] | 100),
blocksize(itime_max/nblock),
iblock(0),
direction(nblock==1? 0:1),
sweeps(0),
sign(1.), 
gf(K_, beta, nblock, blocksize)
//eng_(parms["SEED"] |42), 
//itime_rng(eng_, boost::uniform_int<itime_type>(0,std::numeric_limits<itime_type>::max())), 
//bond_rng(eng_, boost::uniform_int<site_type>(0,n_bond))
//ratio_rng(eng_, boost::uniform_real<>)
{
    //initialize ALPS observables
    initialize_observables();

   if(node==0) {
       print(std::cout); // print parameters to screen 
   }
}


void InteractionExpansion::update()// sweep in one block 
{

   for (unsigned i=0; i< nblock; ++i){

      sweeps++; // one sweep means try go through a block (with steps_per_block updates)
      interaction_expansion_step();                
 
      iblock += direction; 
      //we jump to a new block and calculate gf at its time origin
      //gf.blockjump(iblock*blocksize, tlist, vlist); //from scratch 
      //gf.wrap(iblock*blocksize, tlist, vlist); //fast wrap 
 
      //if hit the end, revert the sweep direction 
      if (iblock == nblock-1 || iblock == 0)
          direction *= -1; 
 
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
