#include "interaction_expansion.hpp"
#include <algorithm>
#include <boost/function.hpp>
#include <boost/bind.hpp>

double npickk( unsigned int n, unsigned int k );  

void InteractionExpansion::Z_to_W2()
{

  if (wormlength!=0) return; 

  std::vector<unsigned int> sites; 
  
  //first site is random 
  //others randomly choose in its neighbor list 
  site_t s = randomint(n_site); 
  sites.push_back(s); 
  unsigned int j= randomint(neighbors[s].size());  
  sites.push_back(neighbors[s][j]);//noway == s because neighbors does not contains the selfsite 
  itime_t tau = beta*random();

  // true means compute_only_weight
  double metropolis_weight = (W2toZ/ZtoW2)* eta2* coef2 *wormcreate_impl(tau, sites, true); // has to be positive 

  /*
  if (metropolis_weight <0.){
     std::cerr << "weight" << metropolis_weight << " "
                                                << std::endl; 

     std::cerr << "operators:" << std::endl ;
     for (int i=0; i< M.creators().size(); ++i){
         std::cerr << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
     }
     std::cerr << std::endl;
  }
  */

  assert(metropolis_weight >-1e-10); // some times we have a problem here 

  if(fabs(metropolis_weight) > random()){
    measurements["ZtoW2"]<< 1.;

    wormcreate_impl(tau, sites, false);

#ifdef VERBOSE 
    std::cout << "create worm at " << "(" << sites[0] << "," << tau << ") "
                                   << "(" << sites[1] << "," << tau << ")" << std::endl; 

    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif

    sign*=metropolis_weight<0.?-1.:1.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); 
    assert(wormlength==2); 

  }else{

    measurements["ZtoW2"]<< 0.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 
    assert(wormlength==0); 
  }
}

void InteractionExpansion::W2_to_Z()
{

    if (wormlength!=2) return; 

    //(1)check whether the worms are within neighbors 
    site_t s = M.creators()[M.num_vertices()*2].s();
    site_t s1 = M.creators()[M.num_vertices()*2+1].s(); 
    if (std::find(neighbors[s].begin(), neighbors[s].end(), s1) == neighbors[s].end())
        return; 

    double metropolis_weight = ZtoW2/(W2toZ*eta2*coef2)* wormdestroy_impl(true, 2);

    assert(metropolis_weight >=0.); 

    if(fabs(metropolis_weight) > random()){ //do the actual update
      measurements["W2toZ"]<< 1.;

      wormdestroy_impl(false, 2);  // false means really perform, not only compute weight

#ifdef VERBOSE 
    std::cout<< "destroy worm" << std::endl; 
    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif

      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()); 
      assert(wormlength==0); 

      sign*=metropolis_weight<0.?-1.:1.;

    }else{
      measurements["W2toZ"]<<0.;

      //do nothing
      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()+wormlength); //that additional 2 is worm 
      assert(wormlength==2); 
    }
}


void InteractionExpansion::W2_to_W4()
{

  if (wormlength!=2) return; 

  std::vector<unsigned int> sites; 
  
  //first site is random 
  //others randomly choose in its neighbor list 
  site_t s = randomint(n_site); 
  sites.push_back(s); 
  unsigned int j= randomint(neighbors[s].size());  
  sites.push_back(neighbors[s][j]);//noway == s because neighbors does not contains the selfsite 
 
  double tau = M.creators()[M.num_vertices()*2].t(); //same tau as exsisting worm 

  // true means compute_only_weight
  double metropolis_weight = (W4toW2/W2toW4)* (coef2*eta4/eta2)* wormcreate_impl(tau, sites, true); // has to be positive 

  /*
  if (metropolis_weight <0.){
     std::cerr << "weight" << metropolis_weight << " "
                                                << std::endl; 

     std::cerr << "operators:" << std::endl ;
     for (int i=0; i< M.creators().size(); ++i){
         std::cerr << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
     }
     std::cerr << std::endl;
  }
  */

  assert(metropolis_weight >-1e-10); // some times we have a problem here 

  if(fabs(metropolis_weight) > random()){
    measurements["W2toW4"]<< 1.;

    wormcreate_impl(tau, sites, false);

#ifdef VERBOSE 
    std::cout << "create worm at " << "(" << sites[0] << "," << tau << ") "
                                   << "(" << sites[1] << "," << tau << ")" << std::endl; 

    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif 
 
    sign*=metropolis_weight<0.?-1.:1.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); 
    assert(wormlength==4); 

  }else{

    measurements["W2toW4"]<< 0.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+2); 
    assert(wormlength==2); 
  }
}

void InteractionExpansion::W4_to_W2()
{
    if (wormlength!=4) return; 

    //remove worm tail 
    //check whether the worms are within neighbors 
    site_t s = M.creators()[M.num_vertices()*2+2].s();
    site_t s1 = M.creators()[M.num_vertices()*2+3].s(); 
    if (std::find(neighbors[s].begin(), neighbors[s].end(), s1) == neighbors[s].end())
        return; 
    
    double metropolis_weight = (W2toW4/W4toW2) * eta2/(eta4*coef2)* wormdestroy_impl(true, 2);

    assert(metropolis_weight >=0.); 

    if(fabs(metropolis_weight) > random()){ //do the actual update
      measurements["W4toW2"]<< 1.;

      wormdestroy_impl(false, 2);  // false means really perform, not only compute weight
   
#ifdef VERBOSE 
    std::cout<< "destroy worm" << std::endl; 
    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif

      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()+2); 
      assert(wormlength==2); 

      sign*=metropolis_weight<0.?-1.:1.;

    }else{
      measurements["W4toW2"]<<0.;

      //do nothing
      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()+wormlength); //that additional 2 is worm 
      assert(wormlength==4); 
    }
}

void InteractionExpansion::Z_to_W4()
{

  if (wormlength!=0) return; 

  std::vector<unsigned int> sites; 
  
  //first site is random 
  //others randomly choose in its neighbor list 
  site_t s = randomint(n_site); 
  sites.push_back(s); 
  for (unsigned int i=0; i<3; ++i){
     unsigned int j= randomint(neighbors[s].size());
     sites.push_back(neighbors[s][j]);//noway == s because neighbors does not contains the selfsite
  }

  itime_t tau = beta*random();

  // true means compute_only_weight
  double metropolis_weight = (W4toZ/ZtoW4)* eta4* coef4 *wormcreate_impl(tau, sites, true); // has to be positive 

  /*
  if (metropolis_weight <0.){
     std::cerr << "weight" << metropolis_weight << " "
                                                << std::endl; 

     std::cerr << "operators:" << std::endl ;
     for (int i=0; i< M.creators().size(); ++i){
         std::cerr << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
     }
     std::cerr << std::endl;
  }
  */

  assert(metropolis_weight >-1e-10); // some times we have a problem here 

  if(fabs(metropolis_weight) > random()){
    measurements["ZtoW4"]<< 1.;

    wormcreate_impl(tau, sites, false);

#ifdef VERBOSE 
    std::cout << "create worm at " << "(" << sites[0] << "," << tau << ") "
                                   << "(" << sites[1] << "," << tau << ")" << std::endl; 

    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif

    sign*=metropolis_weight<0.?-1.:1.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); 
    assert(wormlength==4); 

  }else{

    measurements["ZtoW4"]<< 0.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()); 
    assert(wormlength==0); 
  }
}

void InteractionExpansion::W4_to_Z()
{

    if (wormlength!=4) return; 

    //(1)check whether the worms are within neighbors 
    site_t s = M.creators()[M.num_vertices()*2].s();
    for (unsigned int i=0; i<3; ++i){
        site_t s1 = M.creators()[M.num_vertices()*2+i+1].s();
        if (std::find(neighbors[s].begin(), neighbors[s].end(), s1) == neighbors[s].end())
            return;
    }

    double metropolis_weight = (ZtoW4/W4toZ)/(eta4*coef4)* wormdestroy_impl(true, 4);

    assert(metropolis_weight >=0.); 

    if(fabs(metropolis_weight) > random()){ //do the actual update
      measurements["W4toZ"]<< 1.;

      wormdestroy_impl(false, 4);  // false means really perform, not only compute weight

#ifdef VERBOSE 
    std::cout<< "destroy worm" << std::endl; 
    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif

      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()); 
      assert(wormlength==0); 

      sign*=metropolis_weight<0.?-1.:1.;

    }else{
      measurements["W4toZ"]<<0.;

      //do nothing
      assert(M.creators().size() == M.matrix().rows()); 
      assert(M.creators().size() == 2*M.num_vertices()+wormlength); //that additional 2 is worm 
      assert(wormlength==4); 
    }
}


void InteractionExpansion::wormadd()
{
    if (wormlength==0) return; 

    unsigned int pert_order = M.num_vertices(); 
    //we want to add n vertices  
    unsigned int n = randomint(n_max)+1;
  
    if(pert_order+n > max_order) 
      return; 
  
    std::vector<site_t> sites; 
    std::vector<double> taus; 
  
    itime_t tau; 
    alps::graph_helper<>::bond_descriptor b; 
   
    for (unsigned int i=0; i< n; ++i){
      tau = beta*random(); 
      taus.push_back(tau); 
      taus.push_back(tau); 
  
      b = lattice.bond(randomint(n_bond));
      sites.push_back(lattice.source(b));
      sites.push_back(lattice.target(b));
    }
  
    // true means compute_only_weight
    double metropolis_weight = pow(-beta*V*n_bond, n)/npickk(M.num_vertices()+n, n)*wormadd_impl(taus, sites, true);

  assert(metropolis_weight >=0.); 

  if(fabs(metropolis_weight) > random()){

    std::stringstream obs_name;
    obs_name<<"WormAdd_"<<n;
    measurements[obs_name.str().c_str()] << 1.;

    wormadd_impl(taus, sites, false);

#ifdef VERBOSE
    std::cout << "insert vertex at " << "(" << sites[0] << "," << tau << ") "
                                     << "(" << sites[1] << "," << tau << ")" << std::endl; 
    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif 
 
    sign*=metropolis_weight<0.?-1.:1.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); 
    assert(wormlength>0); 

  }else{

    std::stringstream obs_name;
    obs_name<<"WormAdd_"<<n;
    measurements[obs_name.str().c_str()] << 0.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); 
    assert(wormlength>0); 
  }
}

void InteractionExpansion::wormremove()
{//remove vertex in the presence of worm  

    if (wormlength==0) return; 

    unsigned int pert_order = M.num_vertices(); 
    //we want to remove n vertices  
    unsigned int n = randomint(n_max)+1;

    if(pert_order < n) // not enough vertex to remove (although we still have a worm )
      return;    
    
   //pick up n different vertices  
    std::vector<unsigned int> vertices;
    while (vertices.size()<n){
        unsigned int vertex_nr=randomint(pert_order);// pickup a random vertex 
        if (std::find(vertices.begin(), vertices.end(), vertex_nr) == vertices.end())
            vertices.push_back(vertex_nr); 
    }
    std::sort (vertices.begin(), vertices.end());

    double metropolis_weight = npickk(pert_order, n)/pow(-beta*V*n_bond, n) * wormremove_impl(vertices, true);


    /*
    if (metropolis_weight <-1e-10){
       std::cerr << "weight, pert_order, vertex_nr, lambda:" << metropolis_weight << " "
                                                             << pert_order << " "
                                                             << vertex_nr << " "
                                                             << lambda << " "
                                                             << std::endl; 

       std::cerr << "operators:" << std::endl ;
       for (int i=0; i< M.creators().size(); ++i){
           std::cerr << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
       }
       std::cerr << std::endl;

       abort();  
    }
    */

  assert(metropolis_weight >-1e-10); // some times we have a problem here 

  if(fabs(metropolis_weight)> random()){

    std::stringstream obs_name;
    obs_name<<"WormRemoval_"<<n;
    measurements[obs_name.str().c_str()] << 1.;

    wormremove_impl(vertices, false);

#ifdef VERBOSE 
    std::cout<< "remove vertex:" << std::endl ;
    std::ostream_iterator<unsigned int> out_it (std::cout ,", ");
    std::copy (vertices.begin(), vertices.end(), out_it);

    for (int i=0; i< M.creators().size(); ++i){
        std::cout<< "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout<< std::endl;
#endif
 
    sign*=metropolis_weight<0.?-1.:1.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); //that additional 2 is worm 
    assert(wormlength>0); 

  }else{

    std::stringstream obs_name;
    obs_name<<"WormRemoval_"<<n;
    measurements[obs_name.str().c_str()] << 0.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); 
    assert(wormlength>0); 
  }
}

void InteractionExpansion::wormshift()
{
    if (wormlength==0) return; 


    //(1)
    //select new worm position randomly 
    //std::vector<unsigned int> sites; 
    //for (unsigned int i=0; i<Morder; ++i){
    //    sites.push_back((unsigned int)(random()*n_site));
    //}
    //itime_t tau = beta*random();
    
    //(2)
    //pick up a random site in the worm, then shift it to its neighbor 
    std::vector<unsigned int> sites; 
    for (unsigned int i=0; i<wormlength; ++i){
       sites.push_back(M.creators()[M.num_vertices()*2+i].s());
    }
    unsigned int i= randomint(wormlength);
    site_t s = M.creators()[M.num_vertices()*2+i].s();

    //select a random site in its neighbor list 
    unsigned int j= randomint(neighbors[s].size());  
    sites[i] = neighbors[s][j];

    //random shuffle in the worm, this helps when reduce from M4 to M2 
    boost::function< unsigned int(const unsigned int)> func = boost::bind(&InteractionExpansion::randomint, this, _1); 
    std::random_shuffle (sites.begin(), sites.end(), func); 

    //all sites are identical 
    //if we do this, not to forget to add 1/Ns to the finial answer of M2
    //bool allequal = 
    //   std::find_if(sites.begin() + 1, sites.end(), std::bind1st(std::not_equal_to<site_t>(), sites.front())) == sites.end();
    //if (allequal) return;  
    
    /* 
    //(3)only for Morder = 2 
    //pick up a random site in the worm, then shift it to its neighbor shell 
    assert(Morder ==2); 
    unsigned int i= (unsigned int) (random()*Morder);// 0 or 1 
    site_t si = M.creators()[M.num_vertices()*2+i].s();       //head 
    site_t sj = M.creators()[M.num_vertices()*2+(1-i)].s();   //tail 

    unsigned int dist = disttable(si,sj); 

    DistanceMap dmap = distmap[si]; 
    std::vector<site_t> shell; // the candidates sites with dist-1 and dist+1
    if (dmap.find(dist+1) != dmap.end())
        std::copy (dmap[dist+1].begin(), dmap[dist+1].end(), std::back_inserter(shell)); 

    if (dmap.find(dist-1) != dmap.end())
        std::copy (dmap[dist-1].begin(), dmap[dist-1].end(), std::back_inserter(shell)); 

    assert(shell.size() == neighborshellsize[dist]); 
    sj = shell[(int)(random()*neighborshellsize[dist])]; //shift tail

    unsigned int distnew = disttable(si,sj);

    //std::cout << "shift from " << dist << " to " << distnew << std::endl; 
    
    //accout for imblance of forward/backward moves 
    double propratio = neighborshellsize[dist]/(double)neighborshellsize[distnew]; 

    std::vector<site_t> sites; 
    sites.push_back(si); 
    sites.push_back(sj); 
    */

    //shift time 
    //itime_t tau = beta*random(); //completely random new tau 
    itime_t tau = M.creators()[M.num_vertices()*2].t() + 0.1*beta*(random()-0.5);
    if (tau>beta) tau -= beta; 
    if (tau<0.) tau += beta; 

    // true means compute_only_weight
   double metropolis_weight = wormshift_impl(tau, sites, true); // has to be positive 
   assert(metropolis_weight >-1e-10);  

  if(fabs(metropolis_weight) > random()){
    measurements["WormShift"]<< 1.;
    wormshift_impl(tau, sites, false);

#ifdef VERBOSE  
    std::cout << "shift worm to " << "(" << sites[0] << "," << tau << ") "
                                  << "(" << sites[1] << "," << tau << ")" << std::endl; 
    std::cout << "operators:" << std::endl ;
    for (int i=0; i< M.creators().size(); ++i){
       std::cout << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
    }
    std::cout << std::endl; 
#endif 

    sign*=metropolis_weight<0.?-1.:1.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength); 
    assert(wormlength>0); 

  }else{

    measurements["WormShift"]<<0.;

    assert(M.creators().size() == M.matrix().rows()); 
    assert(M.creators().size() == 2*M.num_vertices()+wormlength);  
    assert(wormlength>0); 
  }
}


