#include "green_function.h"
#include <iostream>

int main(){
    
    double beta = 4.; 
    unsigned nsite = 8; 
    Mat K = Mat::Zero(nsite, nsite); 

    for (unsigned i=0; i< nsite; ++i){
        K(i, (i+1)%nsite) = -1.0; 
        K((i+1)%nsite, i) = -1.0; 
    }

    Green_function gf(nsite, K, beta); 

    double tau1 = 2.513; 
    double tau2 = 0.513;  

    tlist_type tlist; 
    vlist_type vlist; 
    
    for (double tau = 0.1; tau < 4.0; tau += 0.1) {
        tlist.push_back(tau); 
        vlist[tau] = std::make_pair(1, 2); 
    }

    //Mat B = gf.B(tau1, tau2, tlist, vlist); 
    //std::cout << "B:\n" << B << std::endl; 

    Mat G = gf.G(0.,tlist, vlist); 
    std::cout << "G:\n" << G << std::endl; 

    std::cout << "w: " << G.determinant() << std::endl; 

    return 0; 
}
