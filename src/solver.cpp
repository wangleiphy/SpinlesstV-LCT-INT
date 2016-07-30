#include "interaction_expansion.hpp"

///The basic updates: adding and removing vertices.
void InteractionExpansion::interaction_expansion_step()
{
       
   for (unsigned step =0; step < steps_per_block; ++step){

       double update_type=random();

       if(update_type < probs[0]){     
           //std::cout << "before add" << std::endl; 
           add();
           //std::cout << "after add" << std::endl; 
       }else if(update_type < probs[1]){
           //std::cout << "before remove" << std::endl; 
           remove(); 
           //std::cout << "after remove" << std::endl; 
       }else{
           //std::cout << "before shift" << std::endl; 
           shift(); 
           //std::cout << "after shift" << std::endl; 
       }
   }
}
