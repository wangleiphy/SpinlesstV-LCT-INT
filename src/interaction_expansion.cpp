#include "interaction_expansion.hpp"
#include <ctime>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <boost/lexical_cast.hpp>
#include "buildK.h"
#include "bgl.hpp"

InteractionExpansion::InteractionExpansion(alps::params &parms, int node)
:alps::mcbase(parms,node),
lattice(make_deprecated_parameters(parms)),
max_order(boost::lexical_cast<unsigned int>(parms["MAX_ORDER"])),
n_site(lattice.num_sites()),
n_bond(lattice.num_bonds()),
mc_steps((boost::uint64_t)parms["SWEEPS"]),
therm_steps((unsigned long)parms["THERMALIZATION"]),        
beta(boost::lexical_cast<double>(parms["BETA"])), //total projection time 
V(boost::lexical_cast<double>(parms["V"])),                        
tlist(), 
vlist(), 
recalc_period(parms["RECALC_PERIOD"]),
measurement_period(parms["MEASUREMENT_PERIOD"]),
itime_max(boost::lexical_cast<itime_type>(parms["ITIME_MAX"])), 
nblock(boost::lexical_cast<itime_type>(parms["NBLOCKS"])),
steps_per_block(parms["STEPS_PER_BLOCK"] | 100),
blocksize(itime_max/nblock),
iblock(0),
direction(nblock==1? 0:1),
cycles(0), 
sweeps(0),
sign(1.), 
timestep(beta/itime_max), 
window_tau(boost::lexical_cast<time_type>(parms["WINDOWSIZE"])), 
window_upper(static_cast<itime_type>((0.5*beta+0.5*window_tau)/timestep)),//the upper and lower indices  
window_lower(static_cast<itime_type>((0.5*beta-0.5*window_tau)/timestep)),//for the window in the center 
K_(buildK(lattice, boost::lexical_cast<std::string>(parms["BCmodifier"]))), // K of the true Ham 
Ktrial_(buildKtrial(lattice)), // to generate the trial wave function  
gf(lattice, parms["L"], parms["W"], K_, Ktrial_, timestep, itime_max, nblock, blocksize, parms["WRAP_REFRESH_PERIOD"]), 
distmap(get_distmap(lattice)), 
disttable(get_disttable(distmap, n_site)), 
shellsize(get_shellsize(distmap)), 
Add(boost::lexical_cast<double>(parms["Add"])),
Remove(boost::lexical_cast<double>(parms["Remove"])),
probs(), // empty vector 
MEASURE_M4(parms.defined("MEASURE_M4") && boost::lexical_cast<bool>(parms["MEASURE_M4"]))
{   

   probs.push_back(Add); 
   probs.push_back(Add+Remove); 

   if (Add+Remove> 1.0) {
       std::cerr << "Add + Remove > 1.0" << std::endl; 
       abort(); 
   }
   
   //initialize ALPS observables
   initialize_observables();

   if(node==0) {
       print(std::cout); // print parameters to screen 
   }
}

void InteractionExpansion::initialize_tvlist(){

   unsigned Nv = static_cast<unsigned>(0.13*(beta * n_site * V));   //initial number of vertices 
   //std::cout << "Nv: " << Nv << std::endl; 
   for (unsigned i=0; i< Nv; ++i) {
       itime_type itau = randomint(itime_max); 
       tlist.insert(itau); 
  
       std::vector<site_type> sites; 
       alps::graph_helper<>::bond_descriptor b = lattice.bond(randomint(n_bond));
       sites.push_back(lattice.source(b));
       sites.push_back(lattice.target(b));
       vlist[itau] = sites;
   }
   gf.init_with_vertex(tlist, vlist); 
}

void InteractionExpansion::update()
{
      sweeps++; // one sweep means try go through a block (with steps_per_block updates)
      interaction_expansion_step();                

      iblock += direction; 

      //if hit the end, revert the sweep direction 
      if (iblock == nblock-1 || iblock == 0){
          direction *= -1; 
          cycles++; 
      }

      //we jump to a new block and calculate gf at its time origin
      gf.wrap(iblock*blocksize, tlist, vlist); //this is necessary because otherwise we might jump over it_ some empty block 

 
      if(sweeps % recalc_period ==0)
         gf.rebuild(tlist, vlist);

}


void InteractionExpansion::measure(){
  if (sweeps > therm_steps && sweeps%measurement_period ==0){
     time_type tau = gf.tau(); 
     if (fabs(tau - 0.5*beta) < window_tau/2.0) //only measure when we are in the center
        measure_observables();
   } 
}

double InteractionExpansion::fraction_completed() const {
    return (sweeps < therm_steps ? 0. : (sweeps - therm_steps)/double(nblock)/double(mc_steps));
}

void InteractionExpansion::save(alps::hdf5::archive & ar) const {
    mcbase::save(ar);

    std::string context = ar.get_context();
    ar.set_context("/simulation/realizations/0/clones/0/checkpoint");
    
    //copy tlist and vlist to vectors 
    std::vector<itime_type> vt(tlist.begin(), tlist.end());
    std::vector<site_type> vi, vj; 
    for (vlist_type::const_iterator it=vlist.begin(); it!=vlist.end(); ++it) {
          vi.push_back(it->second[0]);
          vj.push_back(it->second[1]); 
    }
    
    ar["sweeps"] << sweeps;
    ar["cycles"] << cycles;
    ar["vt"] << vt;
    ar["vi"] << vi;
    ar["vj"] << vj;

    ar.set_context(context);

}

void InteractionExpansion::load(alps::hdf5::archive & ar) {
    mcbase::load(ar);

    std::string context = ar.get_context();
    ar.set_context("/simulation/realizations/0/clones/0/checkpoint");
    
    std::vector<itime_type> vt; 
    std::vector<site_type> vi, vj; 
    
    ar["sweeps"] >> sweeps;
    ar["cycles"] >> cycles;
    ar["vt"] >> vt;
    ar["vi"] >> vi;  
    ar["vj"] >> vj;
    ar.set_context(context);

    //copy vectors to tlist and vlist
    for (unsigned i=0; i< vt.size(); ++i){
        itime_type itau = vt[i]; 
        tlist.insert(itau); 

        std::vector<site_type> sites; 
        sites.push_back(vi[i]);  
        sites.push_back(vj[i]);  
        vlist[itau] = sites;
    }

    gf.init_with_vertex(tlist, vlist); 
}
