#include "green_function.h"
#include <iostream>

int main(){
    
    time_type beta = 4.; 
    site_type nsite = 8; 
    Mat K = Mat::Zero(nsite, nsite); 

    for (site_type i=0; i< nsite; ++i){
        K(i, (i+1)%nsite) = -1.0; 
        K((i+1)%nsite, i) = -1.0; 
    }

    Green_function gf(K, beta); 

    time_type tau1 = 2.513; 
    time_type tau2 = 0.513;  

    tlist_type tlist; 
    vlist_type vlist; 
    
    for (time_type tau = 0.1; tau < 4.0; tau += 0.1) {
        tlist.insert(tau); 

        vlist[tau].push_back(1);  
        vlist[tau].push_back(2);  
    }

    //Mat B = gf.B(tau1, tau2, tlist, vlist); 
    //std::cout << "B:\n" << B << std::endl; 

    Mat G = gf.G(0.,tlist, vlist); 
    std::cout << "G:\n" << G << std::endl; 

    std::cout << "w: " << G.determinant() << std::endl; 

    return 0; 
}
