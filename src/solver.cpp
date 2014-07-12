#include "interaction_expansion.hpp"

///The basic updates for the InteractionExpansion algorithm: adding and removing vertices.
///This is the heart of InteractionExpansion's code.
void InteractionExpansion::interaction_expansion_step()
{

    if(random() < 0.5){     
//        std::cout << "before add" << std::endl; 
        add();
//        std::cout << "after add" << std::endl; 
    }else{
//        std::cout << "before remove" << std::endl; 
        remove(); 
//        std::cout << "before M2wormcreate" << std::endl; 
    }
}
