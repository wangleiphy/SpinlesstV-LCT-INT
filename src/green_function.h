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
        Green_function(const Mat& K, const time_type beta)
        :ns_(K.rows())
        ,beta_(beta)
        ,itau_(0)
        ,gtau_()
        {
   
         Eigen::SelfAdjointEigenSolver<Mat> ces;
         ces.compute(K);
   
         wK_ = ces.eigenvalues();  
         uK_ = ces.eigenvectors(); 
        
         //initially gtau_ is noninteracting gf at time itau = 0 
         gtau_ = Mat::Identity(ns_, ns_) + expmK(itime_max);
         gtau_ = gtau_.inverse(); 

        }
        
        //return a reference so update of vertex will afftect it 
        Mat& wrap(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {

            if (itau >= itau_) {
                Mat Bmat = B(itau, itau_, tlist, vlist); 
                itau_ = itau; 
                gtau_ = Bmat * gtau_ * Bmat.inverse();
                return gtau_; 
            }else{
                Mat Bmat = B(itau_, itau, tlist, vlist); 
                itau_ = itau; 
                gtau_ = Bmat.inverse() * gtau_  * Bmat; 
                return gtau_; 
            }
        }

         //equal time Green's function at tau 
         Mat G(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) const {// this is very expansive because of inverse
           Mat res = Mat::Identity(ns_, ns_) + B(itau, 0, tlist, vlist) * B(itime_max, itau, tlist, vlist); 
           return res.inverse(); 
         }


    private:
 
         Mat B(const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist) const { // B(tau1) ... B(tau2)
     
             assert(itau1>=itau2); 
     
             tlist_type::const_iterator lower, upper; 
             lower = std::lower_bound (tlist.begin(), tlist.end(), itau2, std::less_equal<itime_type>()); //equal is exclude 
             upper = std::upper_bound (tlist.begin(), tlist.end(), itau1); 

             //std::cout << (lower- tlist.begin())  << " " <<  (upper- tlist.begin())  << " " << tlist.size()   << std::endl;    
             //std::copy(lower, upper, std::ostream_iterator<double>(std::cout, " "));
             //std::cout << std::endl;  
            
             //std::cout << "tlist: " << std::endl; 
             //std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
             //std::cout << std::endl; 

             //std::cout << "tau1, tau2: "<< itau1 << " " << itau2  << std::endl;    
             //std::copy(lower, upper, std::ostream_iterator<itime_type>(std::cout, " "));
             //std::cout << std::endl; 
             
             //upper > tau1 > lower > tau2 
             if (lower == upper ) {// there is no vertex in between tau1 and tau2 
                 return expmK(itau1 - itau2); 
             }else{

                 Mat res = expmK(*lower - itau2);
                 for (tlist_type::const_iterator it1 =lower, it2 =++lower; it1!=upper; ++it1, ++it2) {
                    
                     itime_type itau = *it1; 
                     res.row(vlist[itau][0]) *= -1.; 
                     res.row(vlist[itau][1]) *= -1.; 
                     //std::cout << "act vertex " << std::endl; 
                 
                     itime_type ditau = (it2 ==upper) ? itau1 - itau: *it2 - itau; 
                 
                     //std::cout << "itau, ditau " << itau << " " << ditau << std::endl; 
                     res = expmK( ditau ) *res; 
                 }
                 //std::cout << "##############" << std::endl; 
                 return res; 
             }
         }
    
     
         Mat expmK(const itime_type itau) const {// exp(-tau * K)
             assert(itau>=0); 

             if (itau ==0)
                 return Mat::Identity(ns_, ns_); 

             Eigen::VectorXd v(ns_); 
             for(unsigned l=0; l<ns_; ++l) 
                 v(l) = exp(-itime2time(itau) * wK_(l)); 
             return uK_ * v.asDiagonal() * uK_.adjoint(); 
         }


        time_type itime2time(const itime_type itau) const{
            return beta_ * itau/itime_max; 
        }

    private:
        const site_type ns_; 
        const time_type beta_; 

        //eigen value and vectors of K 
        Eigen::VectorXd wK_; 
        Mat uK_; 

        itime_type itau_; 
        Mat gtau_; 
};

#endif
