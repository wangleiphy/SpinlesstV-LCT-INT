#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include "types.h"
#include <iterator>
#include <iostream>

class Green_function{

    public:

        ///constructor: how many time slices, how many sites
        Green_function(const Mat& K, const time_type beta, const unsigned nblock, const unsigned blocksize)
        :ns_(K.rows())
        ,beta_(beta)
        ,itau_(0)
        ,U_(Mat::Identity(ns_, ns_))
        ,D_(Mat::Zero(ns_, ns_))
        ,V_(Mat::Identity(ns_, ns_))
        ,nblock_(nblock)// number of blocks 
        ,blocksize_(blocksize)
        {
   
         Eigen::SelfAdjointEigenSolver<Mat> ces;
         ces.compute(K);
   
         wK_ = ces.eigenvalues();  
         uK_ = ces.eigenvectors(); 
         uKdag_ = ces.eigenvectors().adjoint(); 

         //std::cout << "U*Udag:\n" << uK_ * uKdag_ << std::endl; 
        
         //initially gtau_ is noninteracting gf = 1./(1+ exp(-E*beta))
         for(site_type l=0; l<ns_; ++l) 
            D_(l, l) = wK_(l)>0. ? 1./(1.+exp(-beta*wK_(l))) : exp(beta_*wK_(l))/(1.+exp(beta_*wK_(l))) ; // it is samething, to avoid overflow 
        }

        void rebuild(const tlist_type& tlist, vlist_type& vlist){

            Mat gtau = Gstable(itau_, tlist, vlist); // from scratch 

            double max_diff = ((gtau - U_*D_*V_).cwiseAbs()).maxCoeff(); 
            if(max_diff > 1.e-6)
              std::cout<<"WARNING: roundoff errors " <<max_diff << std::endl;
           
            Eigen::JacobiSVD<Mat> svd(gtau, Eigen::ComputeThinU | Eigen::ComputeThinV); 

            U_ = svd.matrixU(); 
            D_ = svd.singularValues().asDiagonal(); 
            V_ = svd.matrixV().adjoint();
        }
        
        /*
        //jump to a new time itau from scratch
        void blockjump(const itime_type itau, const tlist_type& tlist, vlist_type& vlist){
             gtau_ = G(itau, tlist, vlist);  
             itau_ = itau; 
        }
        */

        Mat gtau() const {
            return U_* D_* V_; 
        }

        double gij(const site_type si, const site_type sj)const {
            // (U gtau U^{dagger} )_ij 
            // gtau in eigen basis is stored as U_*D_*V_ 
            return  (uK_.row(si) *U_) *D_ * (V_ * uKdag_.col(sj));  
        }

        Eigen::VectorXd denmat(const site_type si)const {
            // (U gtau U^{dagger} )_ij 
            // gtau in eigen basis is stored as U*D*V 
            return  uK_ *(U_ *D_ * (V_ * uKdag_.col(si)));  
        }
            
        //update changes D_ and V_ 
        void update(const site_type si, const site_type sj, const double gij){

             //update D_
             D_ -=  ( (D_* (V_* uKdag_.col(sj))) * ((uK_.row(si)*U_)*D_)
                    + (D_* (V_* uKdag_.col(si))) * ((uK_.row(sj)*U_)*D_)
                    )/gij; 
             
             // * U^\dagger V U from right 
             V_ -= 2.* (V_* uKdag_.col(si)) * uK_.row(si) + 2.* (V_*uKdag_.col(sj)) * uK_.row(sj); // this thing is fixed 

             Eigen::JacobiSVD<Mat> svd(D_, Eigen::ComputeThinU | Eigen::ComputeThinV); 
             
             U_ = U_*svd.matrixU(); 
             D_ = svd.singularValues().asDiagonal(); 
             V_ = svd.matrixV().adjoint()*V_;
        }
        
        //wrap does not change D_  
        void wrap(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {

            if (itau >= itau_) {
                // B G B^{-1}
                propagator1(-1, itau, itau_, tlist, vlist, U_);  // B(tau1) ... B(tau2) *U_  
                propagator1(1, itau, itau_, tlist, vlist, V_); // V_ * B^{-1}(tau2) ... B^{-1}(tau1)

                itau_ = itau; 

            }else{

                // B^{-1} G B 
                propagator2(1, itau_, itau, tlist, vlist, U_); //  B^{-1}(tau2) ... B^{-1}(tau1) * U_
                propagator2(-1, itau_, itau, tlist, vlist, V_);   //  V_ * B(tau1) ... B(tau2)

                itau_ = itau; 
            }
        }

         //equal time Green's function at tau 
         Mat G(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {// this is very expansive because of inverse
           //Mat res = Mat::Identity(ns_, ns_) + B(itau, 0, tlist, vlist) * B(itime_max, itau, tlist, vlist); 
           //Mat res = Mat::Identity(ns_, ns_) + B_tau_0 * B_beta_tau; 
       
            Mat B_tau_0 = Mat::Identity(ns_, ns_); 
            propagator1(-1, itau, 0, tlist, vlist, B_tau_0);
            
            Mat B_beta_tau =  Mat::Identity(ns_, ns_);
            propagator2(-1, itime_max, itau, tlist, vlist, B_beta_tau);
            
            Mat res = Mat::Identity(ns_, ns_) + B_tau_0 * B_beta_tau; 
            return res.inverse(); 
         }

        Mat Gstable(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {// this is very expansive because of inverse
          //B_tau_0 = U1*D1*V1
          Mat U1 = Mat::Identity(ns_, ns_); 
          Mat D1 = Mat::Identity(ns_, ns_); 
          Mat V1 = Mat::Identity(ns_, ns_); 
           
          unsigned b = itau/blocksize_; //block index 

           for (unsigned ib=0; ib< b; ++ib) {
                propagator1(-1, (ib+1)*blocksize_, ib*blocksize_, tlist, vlist, U1);

                Eigen::JacobiSVD<Mat> svd(U1*D1, Eigen::ComputeThinU | Eigen::ComputeThinV); 
                U1 = svd.matrixU();
                D1 = svd.singularValues().asDiagonal();
                V1 = svd.matrixV().adjoint()* V1;
           }
           propagator1(-1, itau, b*blocksize_, tlist, vlist, U1);
        
          //B_beta_tau = U2*D2*V2
          Mat U2 = Mat::Identity(ns_, ns_); 
          Mat D2 = Mat::Identity(ns_, ns_); 
          Mat V2 = Mat::Identity(ns_, ns_); 

          for (unsigned ib=nblock_; ib> b+1; --ib) {
            propagator2(-1, ib*blocksize_, (ib-1)*blocksize_, tlist, vlist, V2);

            Eigen::JacobiSVD<Mat> svd(D2*V2, Eigen::ComputeThinU | Eigen::ComputeThinV); 
            U2 = U2*svd.matrixU();
            D2 = svd.singularValues().asDiagonal();
            V2 = svd.matrixV().adjoint();
          }

           propagator2(-1, (b+1)*blocksize_, itau, tlist, vlist, V2);

           Mat res= U1.inverse()*V2.inverse() + D1 * V1 * U2 * D2;

           Eigen::JacobiSVD<Mat> svd(res, Eigen::ComputeThinU | Eigen::ComputeThinV); 
           
           Mat U, D, V;   
           U = svd.matrixU(); 
           D = svd.singularValues().asDiagonal(); 
           V = svd.matrixV().adjoint();

           return (V*V2).inverse() * D.inverse() * (U1*U).inverse(); 
         }
        
         /*
         //U and U^{dgger} 
         const Mat& U() const{
             return uK_;  
         }

         const Mat& Udag() const{
             return uKdag_;  
         }
         */

        
         // it can do  B(tau1)... B(tau2) * A  when sign = -1
         // or        A* B(tau2)^{-1} ... B(tau1)^{-1} when sign = 1 
         void propagator1(const int sign, const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist, Mat& A) { 
     
             assert(itau1>=itau2); 

             std::string side = (sign ==-1) ? "L" : "R"; 
     
             tlist_type::const_iterator lower, upper; 
             lower = std::lower_bound (tlist.begin(), tlist.end(), itau2, std::less_equal<itime_type>()); //equal is exclude 
             upper = std::upper_bound (tlist.begin(), tlist.end(), itau1); 

             //std::cout << "vertices at" << std::endl; 
             //std::copy(lower, upper, std::ostream_iterator<itime_type>(std::cout, " "));
             //std::cout << std::endl;  

             //upper > tau1 > lower > tau2 
             if (lower == upper ) {// there is no vertex in between tau1 and tau2 
                 return Kprop(sign, itau1 - itau2, side, A); 

             }else{

                 Kprop(sign, *lower - itau2, side, A);

                 //std::cout << "initial " << *lower << " " << itau2 << " " <<  *lower - itau2 << std::endl; 
                 for (tlist_type::const_iterator it1 =lower, it2 =++lower; it1!=upper; ++it1, ++it2) {
                    
                     itime_type itau = *it1; 
                     Vprop(vlist[itau][0], vlist[itau][1], side, A); 
                     //std::cout << "act vertex at " << itau  << std::endl; 
                 
                     itime_type ditau = (it2 ==upper) ? itau1 - itau: *it2 - itau; 
                     Kprop(sign, ditau , side, A); 
                     //std::cout << "Kprop " <<   ((it2 ==upper) ? itau1: *it2)  << " " << itau  << " " << ditau << std::endl; 
                 }
                 //std::cout << "##############" << std::endl; 
             }
         }

         // it can do A* B(tau1)... B(tau2)  when sign = -1
         // or        B(tau2)^{-1} ... B(tau1)^{-1}* A when sign = 1 
         void propagator2(const int sign, const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist, Mat& A) { 
     
             assert(itau1>=itau2); 

             std::string side = (sign ==1) ? "L" : "R"; 
     
             tlist_type::const_iterator lower, upper; 
             lower = std::lower_bound (tlist.begin(), tlist.end(), itau2, std::less_equal<itime_type>()); 
             upper = std::upper_bound (tlist.begin(), tlist.end(), itau1); 
                
             //std::cout << "vertices at" << std::endl; 
             //std::copy(lower, upper, std::ostream_iterator<itime_type>(std::cout, " "));
             //std::cout << std::endl;  

             //upper > tau1 > lower > tau2 
             if (lower == upper ) {// there is no vertex in between tau1 and tau2 
                 return Kprop(sign, itau1 - itau2, side, A); 

             }else{

                 Kprop(sign, itau1 - *(--upper), side, A);
                 //std::cout << "initial " << *upper << " " << itau1 << " " << itau1 - *upper << std::endl; 
                 for (tlist_type::const_iterator it1 =upper, it2 =--upper; it1!=lower; --it1, --it2) {
                    
                     itime_type itau = *it1; 
                     Vprop(vlist[itau][0], vlist[itau][1], side, A); 
                     //std::cout << "act vertex at " << itau  << std::endl; 
                 
                     itime_type ditau = itau - *it2; 
                     Kprop(sign, ditau , side, A); 
                     //std::cout << "Kprop " <<  itau << " " << *it2  << " " << ditau << std::endl; 
                 }
                 {
                     //the last step by hand (it1 = lower, it2 point to tau2) 
                     itime_type itau = *lower; 
                     Vprop(vlist[itau][0], vlist[itau][1], side, A); 
                     //std::cout << "act vertex at " << itau  << std::endl; 
             
                     itime_type ditau = itau - itau2; 
                     Kprop(sign, ditau , side, A); 
                     //std::cout << "Kprop " <<  itau << " " <<  itau2 << " " << ditau << std::endl; 
                 }
                 //std::cout << "##############" << std::endl; 
             }
         }
    
         /*
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
         */

    private:

        void Kprop(const int sign, const itime_type itau, const std::string side, Mat& A) const{
            if (itau ==0) 
                return; 
                
            time_type tau = itime2time(itau); 
            //std::cout << "tau: " << tau << std::endl; 

            if (side == "L") { //  exp( sign * tau *w) * A 
                for(unsigned l=0; l<ns_; ++l) 
                  A.row(l) *= exp( sign * tau * wK_(l));
            }else if (side == "R") { // A * exp(sign * tau*w)
                for(unsigned l=0; l<ns_; ++l) 
                  A.col(l) *= exp( sign * tau * wK_(l));
            }
        }

        void Vprop(const site_type si, const site_type sj, const std::string side,  Mat& A) const{ // A*U^{dagger} V U or U^{dagger} V U * A 
            
            if (side == "R"){
                A -= 2.* (A*uKdag_.col(si))* uK_.row(si) + 2.* (A*uKdag_.col(sj)) * uK_.row(sj);
            }else if (side == "L"){
                A -= 2.* uKdag_.col(si) * (uK_.row(si)* A) + 2.* uKdag_.col(sj) * (uK_.row(sj)* A); 
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

        //gtau = U_*D_*V_
        Mat U_, D_, V_; 
        
        //blocks are used when calculating G from scratch 
        unsigned nblock_, blocksize_; 

};

#endif
