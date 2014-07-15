#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include "types.h"
#include <iterator>
#include <iostream>

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
         uKdag_ = ces.eigenvectors().adjoint(); 

         //std::cout << "U*Udag:\n" << uK_ * uKdag_ << std::endl; 
        
         //initially gtau_ is noninteracting gf at time itau = 0 
         gtau_ = Mat::Identity(ns_, ns_) + expmK(-1, itime_max); // 1+ exp(-K*beta)
         gtau_ = gtau_.inverse(); 

        }

        void rebuild(const tlist_type& tlist, vlist_type& vlist){

            Mat gtau = G(itau_, tlist, vlist); 

            double max_diff = ((gtau - gtau_).cwiseAbs()).maxCoeff(); 
            if(max_diff > 1.e-8)
              std::cout<<"WARNING: roundoff errors " <<max_diff << std::endl;
            
            gtau_ = gtau; 
        }
        
        //jump to a new time itau from scratch
        void blockjump(const itime_type itau, const tlist_type& tlist, vlist_type& vlist){
             gtau_ = G(itau, tlist, vlist);  
             itau_ = itau; 
        }

        const Mat& gtau()const {return gtau_; }
        
        //return a reference so update of vertex will afftect it 
        Mat& wrap(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {


            if (itau >= itau_) {

                //std::cout << "B*Binv: " << itau << " " << itau_ << "\n" << B(itau, itau_, tlist, vlist) * Binv(itau, itau_, tlist, vlist) << std::endl;

                //std::cout << "det(B):\n" <<  B(itau, itau_, tlist, vlist).determinant() << std::endl; 

                //std::cout << "Binv:\n" <<  Binv(itau, itau_, tlist, vlist) << std::endl; 
                //std::cout << "B^{-1}:\n" <<  B(itau, itau_, tlist, vlist).inverse() << std::endl; 

                gtau_ = B(itau, itau_, tlist, vlist) * gtau_ * Binv(itau, itau_, tlist, vlist);
                itau_ = itau; 
            }else{

                //std::cout << "B*Binv: " << itau << " " << itau_ << "\n" << B(itau_, itau, tlist, vlist) * Binv(itau_, itau, tlist, vlist) << std::endl;

                //std::cout << "det(B):\n" <<  B(itau_, itau, tlist, vlist).determinant() << std::endl; 
                //std::cout << "Binv:\n" <<  Binv(itau_, itau, tlist, vlist) << std::endl; 
                //std::cout << "B^{-1}:\n" <<  B(itau_, itau, tlist, vlist).inverse() << std::endl; 

                gtau_ = Binv(itau_, itau, tlist, vlist) * gtau_  * B(itau_, itau, tlist, vlist); 
                itau_ = itau; 
            }
                return gtau_; 
        }

         //equal time Green's function at tau 
         Mat G(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) const {// this is very expansive because of inverse
           Mat res = Mat::Identity(ns_, ns_) + B(itau, 0, tlist, vlist) * B(itime_max, itau, tlist, vlist); 
           return res.inverse(); 
         }
        
         //U and U^{dgger} 
         const Mat& U() const{
             return uK_;  
         }

         const Mat& Udag() const{
             return uKdag_;  
         }

         Mat B(const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist) const { // B(tau1) ... B(tau2)
             return BorBinv(-1, itau1, itau2, tlist, vlist); 
         }

         Mat Binv(const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist) const { // B^{-1}(tau2) ... B^{-1}(tau1)
             return BorBinv(1, itau1, itau2, tlist, vlist); 
         }

    private:
             

         Mat BorBinv(const int sign, const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist) const { 
     
             assert(itau1>=itau2); 
     
             tlist_type::const_iterator lower, upper; 
             lower = std::lower_bound (tlist.begin(), tlist.end(), itau2, std::less_equal<itime_type>()); //equal is exclude 
             upper = std::upper_bound (tlist.begin(), tlist.end(), itau1); 
             
             //upper > tau1 > lower > tau2 
             if (lower == upper ) {// there is no vertex in between tau1 and tau2 
                 return expmK(sign, itau1 - itau2); 
             }else{

                 Mat res = expmK(sign, *lower - itau2);
                 //std::cout << "initial " << *lower << " " << itau2 << " " <<  *lower - itau2 << std::endl; 
 
                 for (tlist_type::const_iterator it1 =lower, it2 =++lower; it1!=upper; ++it1, ++it2) {
                    
                     itime_type itau = *it1; 
                     Vprop(sign, vlist[itau][0], vlist[itau][1], res); 
                     //std::cout << "act vertex at " << itau  << std::endl; 
                 
                     itime_type ditau = (it2 ==upper) ? itau1 - itau: *it2 - itau; 
                     Kprop(sign, ditau , res); 
                     //std::cout << "Kprop " <<   ((it2 ==upper) ? itau1: *it2)  << " " << itau  << " " << ditau << std::endl; 
                 }
                 //std::cout << "##############" << std::endl; 
                 return res; 
             }
         }
    
     
         Mat expmK(const int sign, const itime_type itau) const {// exp(sign * tau * w)
             assert(itau>=0); 

             if (itau ==0)
                 return Mat::Identity(ns_, ns_); 
            
             time_type tau = itime2time(itau); 

             Eigen::VectorXd v(ns_); 
             for(unsigned l=0; l<ns_; ++l) 
                 v(l) = exp(sign * tau * wK_(l)); 
             
             return v.asDiagonal(); 
         }

        void Kprop(const int sign, const itime_type itau, Mat& A) const{// exp(-tau *w) * A  or A * exp(tau*w)
            if (itau ==0) 
                return; 
                
            time_type tau = itime2time(itau); 
            //std::cout << "tau: " << tau << std::endl; 

            if (sign ==-1) {
                for(unsigned l=0; l<ns_; ++l) 
                  A.row(l) *= exp(- tau * wK_(l));
            }else if (sign == 1) {
                for(unsigned l=0; l<ns_; ++l) 
                  A.col(l) *= exp( tau * wK_(l));
            }
        }

        void Vprop(const int sign, const site_type si, const site_type sj , Mat& A) const{ // sign=1: A*U^{dagger} V U;  sign=-1 : U^{dagger} V U * A 
            
            if (sign == 1){
                A -= 2.* (A*uKdag_.col(si))* uK_.row(si) + 2.* (A*uKdag_.col(sj)) * uK_.row(sj);
                //A = A*uKdag_; 
                //A.col(si) *= -1; 
                //A.col(sj) *= -1; 
                //A = A* uK_; 

            }else if (sign== -1){
                A -= 2.* uKdag_.col(si) * (uK_.row(si)* A) + 2.* uKdag_.col(sj) * (uK_.row(sj)* A); 
                //A = uK_*A; 
               //A.row(si) *= -1; 
               // A.row(sj) *= -1; 
                //A = uKdag_ * A; 
            }
        }

        time_type itime2time(const itime_type itau) const{
            return beta_*itau/itime_max; 
        }

    private:
        const site_type ns_; 
        const time_type beta_; 

        //eigen value and vectors of K 
        Eigen::VectorXd wK_; 
        Mat uK_; 
        Mat uKdag_; 

        itime_type itau_; 
        Mat gtau_; 
};

#endif
