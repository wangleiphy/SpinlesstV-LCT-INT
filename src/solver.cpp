#include "interaction_expansion.hpp"

///The basic updates for the InteractionExpansion algorithm: adding and removing vertices.
///This is the heart of InteractionExpansion's code.
void InteractionExpansion::interaction_expansion_step()
{
       
   for (unsigned step =0; step < steps_per_block; ++step){
       if(random() < 0.5){     
           //std::cout << "before add" << std::endl; 
           add();
           //std::cout << "after add" << std::endl; 
       }else{
           //std::cout << "before remove" << std::endl; 
           remove(); 
           //std::cout << "after remove" << std::endl; 
       }
   }
}
