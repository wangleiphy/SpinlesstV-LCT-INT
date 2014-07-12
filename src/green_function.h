#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include "types.h"
#include <iterator>
#include <iostream>

// G(tau2) = B(tau2, tau1) G(tau1) B^{-1}(tau2, tau1)
//void wrapgf(const double tau1, const double tau2){
//}

class Green_function{

    public:

        ///constructor: how many time slices, how many sites
        Green_function(unsigned int nsite, const Mat& K, const double beta)
        :ns_(nsite)
        ,beta_(beta)
        {
   
         Eigen::SelfAdjointEigenSolver<Mat> ces;
         ces.compute(K);
   
         wK_ = ces.eigenvalues();  
         uK_ = ces.eigenvectors(); 
   
        }
 
         Mat B(const double tau1, const double tau2, const tlist_type& tlist, vlist_type& vlist) const { // B(tau1) ... B(tau2)
     
             assert(tau1>=tau2); 
     
             tlist_type::const_iterator lower, upper; 
             lower = std::lower_bound (tlist.begin(), tlist.end(), tau2); 
             upper = std::upper_bound (tlist.begin(), tlist.end(), tau1); 
            
             //if (lower == upper){
             //    return Mat::Identity(nsite, nsite); 
             //}else{

             std::cout << (lower- tlist.begin())  << " " <<  (upper- tlist.begin())  << " " << tlist.size()   << std::endl;    
             std::cout << *lower << " " <<  *upper  << std::endl;    
             //std::copy(lower, upper, std::ostream_iterator<double>(std::cout, " "));
             //std::cout << std::endl;  
             
             //upper > tau1 > lower > tau2 
             if (lower == upper ) {// there is no vertex in between tau1 and tau2 
                 return expmK(tau1 - tau2); 
             }else{
            
                 Mat res = expmK(*lower - tau2);
                 for (tlist_type::const_iterator it=lower; it!=upper; ++it) {
                    
                     double tau = *it; 
                     res.row(vlist[tau].first) *= -1.; 
                     res.row(vlist[tau].second) *= -1.; 
                 
                     double dtau = (it == upper-1) ? tau1 - tau: *(it+1) - tau; 
                 
                     std::cout << "tau, dtau " << tau << " " << dtau << std::endl; 
                     res = expmK( dtau ) *res; 
                 }
                 return res; 
             }
         }
     
         //equal time Green's function at tau 
         Mat G(const double tau, const tlist_type& tlist, vlist_type& vlist) const {// this is very expansive because of inverse
           Mat res = Mat::Identity(ns_, ns_) + B(tau, 0., tlist, vlist) * B(beta_, tau, tlist, vlist); 
           return res.inverse(); 
         }
     
         Mat expmK(const double tau) const {// exp(-tau * K)
             Eigen::VectorXd v(ns_); 
             for(unsigned l=0; l<ns_; ++l) 
                 v(l) = exp(-tau * wK_(l)); 
             return uK_ * v.asDiagonal() * uK_.adjoint(); 
         }

    private:
        const unsigned ns_; 
        const double beta_; 

        //eigen value and vectors of K 
        Eigen::VectorXd wK_; 
        Mat uK_; 
 
};

#endif
