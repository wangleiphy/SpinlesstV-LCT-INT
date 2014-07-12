#include <alps/ngs.hpp>
#include <alps/ngs/scheduler/parseargs.hpp>
#include "interaction_expansion.hpp"

void InteractionExpansion::test(){

   {//this block adds ONE vertex
      std::vector<site_t> sites;  
      std::vector<double> taus; 
      sites.push_back(1); 
      sites.push_back(14); 
      sites.push_back(12); 
      sites.push_back(13); 

      double tau = 1.90143; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      tau = 1.46399; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      add_impl(taus, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << M.num_vertices() << std::endl; 
   }
   
   {//this block adds vertex
      std::vector<site_t> sites;  
      std::vector<double> taus; 
      sites.push_back(3); 
      sites.push_back(4); 

      sites.push_back(8); 
      sites.push_back(7); 

      double tau = 0.234; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      tau = 0.3; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      add_impl(taus, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << M.num_vertices() << std::endl; 
   }

   {//this block adds vertex
      std::vector<site_t> sites;  
      std::vector<double> taus; 
      sites.push_back(4); 
      sites.push_back(5); 

      sites.push_back(12); 
      sites.push_back(13); 

      double tau = 0.124; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      tau = 0.3321; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      add_impl(taus, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add vertex" << std::endl; 
      std::cout << "number of vertices: " << M.num_vertices() << std::endl; 
   }

   {//this block removes vertex 
      std::vector<unsigned int> vertices; // should be sorted 
      vertices.push_back(1); 
      vertices.push_back(4); 

      remove_impl(vertices, false);
      std::cout << "remove vertex" << std::endl; 
      std::cout << "number of vertices: " << M.num_vertices() << std::endl; 
   }


   {//this block adds worm 
      std::vector<site_t> sites;  
      std::vector<double> taus;  
      sites.push_back(3); 
      sites.push_back(2); 

      double tau = 0.14536; 

      wormcreate_impl(tau, sites, false); //actual update M.matrix() and M.num_vertices(), false means not only calculate weight 
      std::cout << "add worm" << std::endl; 
   }

   {//this block add vertex in the presence of worm  
      std::vector<site_t> sites;  
      sites.push_back(2); 
      sites.push_back(9); 

      std::vector<double> taus; 
      double tau = 0.445; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      wormadd_impl(taus, sites, false);
      std::cout << "add vertex in worm space" << std::endl; 
      std::cout << "number of vertices: " << M.num_vertices() << std::endl; 
   }

   {//this block add vertex in the presence of worm  
      std::vector<site_t> sites;  
      sites.push_back(5); 
      sites.push_back(6); 

      std::vector<double> taus; 
      double tau = 0.323; 
      taus.push_back(tau); 
      taus.push_back(tau); 

      wormadd_impl(taus, sites, false);
      std::cout << "add vertex in worm space" << std::endl; 
      std::cout << "number of vertices: " << M.num_vertices() << std::endl; 
   }

   {//this block remove vertex in the presence of worm  
      std::vector<unsigned int> vertices; // should be sorted 
      vertices.push_back(0); 
      vertices.push_back(1); 
      vertices.push_back(2); 
      vertices.push_back(3); 

      wormremove_impl(vertices, false);

      std::cout << "remove vertex in worm space" << std::endl; 
      std::cout << "number of vertices: " << M.num_vertices() << std::endl; 
   }

   {//this block shifts the worm 
      std::vector<site_t> sites;  
      sites.push_back(3); 
      sites.push_back(2); 
      double tau = 0.321; 

      wormshift_impl(tau, sites, false);
      std::cout << "shift worm" << std::endl; 
   }

   {//this block shifts the worm 
      std::vector<site_t> sites;  
      sites.push_back(1); 
      sites.push_back(1); 
      double tau = 0.1323232; 

      wormshift_impl(tau, sites, false);
      std::cout << "shift worm" << std::endl; 
   }

   {//this block destroy worm
      wormdestroy_impl(false, 2);
      std::cout << "remove worm" << std::endl; 
   }

   std::cerr << "cdag:" << std::endl ;
   for (int i=0; i< M.creators().size(); ++i){
       std::cerr << "(" << M.creators()[i].s() << "," << M.creators()[i].t() << ") ";
   }
   std::cerr << std::endl;
   std::cout << "M.matrix() from fast update:\n" << M.matrix() << std::endl; 
   std::cout << "det(M) from fast update= " << M.matrix().determinant() << std::endl; 

}


int main(int argc, char** argv){
   alps::parseargs options(argc, argv); 
   alps::params params(options.input_file);

   InteractionExpansion sim(params, 0); 
   std::cout << "initialization done" << std::endl; 
   sim.test(); 
   sim.build_matrix(); 
   return 0; 
}

