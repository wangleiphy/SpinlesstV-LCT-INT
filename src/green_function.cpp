#include "green_function.h"
#include <iostream>
#include <boost/random.hpp>

int main(){
    
    time_type beta = 4.; 
    site_type nsite = 8; 
    Mat K = Mat::Zero(nsite, nsite); 

    for (site_type i=0; i< nsite; ++i){
        K(i, (i+1)%nsite) = -1.0; 
        K((i+1)%nsite, i) = -1.0; 
    }

    Green_function gf(K, beta); 


    typedef boost::mt19937 engine_type;
    engine_type eng(42);

    //generator for itime index 
    typedef boost::uniform_int<itime_type> dist_type;
    typedef boost::variate_generator<engine_type&, dist_type> rng_type;
    rng_type itime_rng(eng,dist_type(0,std::numeric_limits<itime_type>::max()));

    tlist_type tlist; 
    vlist_type vlist; 
    
    unsigned Nv = 5;   
    itime_type itau; 
    for (unsigned i=0; i< Nv; ++i) {
        itau = itime_rng(); 

        tlist.insert(itau); 

        vlist[itau].push_back(1);  
        vlist[itau].push_back(2);  
    }

    //itime_type itau1 = 3143890026; //(2.5/beta)*std::numeric_limits<itime_type>::max(); 
    //itime_type itau2 = 1608637542; //(0.513/beta)*std::numeric_limits<itime_type>::max();  

    //Mat B = gf.B(itau1, itau2, tlist, vlist); 
    //std::cout << "B:\n" << B << std::endl; 

    Mat G = gf.G(itau, tlist, vlist); 
    std::cout << "G:\n" << G << std::endl; 
    std::cout << "w: " << G.determinant() << std::endl; 

    return 0; 
}

