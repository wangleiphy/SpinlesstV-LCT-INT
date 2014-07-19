#include "green_function.h"
#include <iostream>
#include <boost/random.hpp>

int main(){
    
    time_type beta = 16.; 
    site_type nsite = 8; 
    Mat K = Mat::Zero(nsite, nsite); 

    for (site_type i=0; i< nsite; ++i){
        K(i, (i+1)%nsite) = -1.0; 
        K((i+1)%nsite, i) = -1.0; 
    }

    unsigned nblock = 5; 
    unsigned blocksize = itime_max/nblock; 
    Green_function gf(K, beta, nblock, blocksize, 10, 10); 

    typedef boost::mt19937 engine_type;
    engine_type eng(42);

    //generator for itime index 
    typedef boost::uniform_int<itime_type> dist_type;
    typedef boost::variate_generator<engine_type&, dist_type> rng_type;
    rng_type itime_rng(eng,dist_type(0,itime_max));

    tlist_type tlist; 
    vlist_type vlist; 
    
    unsigned Nv = 1;   
    itime_type itau; 
    for (unsigned i=0; i< Nv; ++i) {
        itau = itime_rng(); 

        tlist.insert(itau); 

        std::vector<site_type> sites; 
        sites.push_back(1);  
        sites.push_back(2);  

        vlist[itau] = sites;
    }


    //Mat B = gf.B(itau1, itau2, tlist, vlist); 
    //std::cout << "B:\n" << B << std::endl; 

    itau = itime_rng(); 

    Mat G = gf.G(itau, tlist, vlist); 
    std::cout << "G:\n" << G << std::endl; 
        
    Mat Gstable = gf.Gstable(itau, tlist, vlist); 
    std::cout << "Gstable:\n" << Gstable << std::endl; 

    /*
    std::cout << "w: " << G.determinant() << std::endl; 

    std::cout << "tlist: "; 
    std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
    std::cout << std::endl; 
          
    itime_type itau1 = 3143890026; //(2.5/beta)*itime_max; 
    itime_type itau2 = 787846414; //(0.513/beta)*itime_max;  

    Mat B = Mat::Identity(nsite, nsite); 
    gf.propagator1(-1, itau1, itau2, tlist, vlist, B);
    std::cout << "B:\n" << B << std::endl; 

    std::cout << "###########" << std::endl; 
    Mat Binv = Mat::Identity(nsite, nsite); 
    gf.propagator1(1, itau1, itau2, tlist, vlist, Binv);

    std::cout << "Binv:\n" << Binv << std::endl; 
    std::cout << "B*Binv:\n" << B*Binv << std::endl; 
    */

    return 0; 
}

