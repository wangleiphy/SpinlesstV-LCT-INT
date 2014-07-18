#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include "types.h"
#include <boost/tuple/tuple.hpp>
#include <iterator>
#include <iostream>

class Green_function{

    public:

        Green_function(const Mat& K, const time_type beta, const time_type timestep, const unsigned nblock, const itime_type blocksize, const unsigned update_refresh_period, const unsigned wrap_refresh_period)
        :ns_(K.rows())
        ,beta_(beta)
        ,timestep_(timestep)
        ,itau_(0)
        ,U_()
        ,D_()
        ,V_()
        ,blocksize_(blocksize)
        ,update_refresh_counter_(0)
        ,update_refresh_period_(update_refresh_period)
        ,wrap_refresh_counter_(0)
        ,wrap_refresh_period_(wrap_refresh_period)
        ,Lstorage_(nblock)
        ,Rstorage_(nblock)
        {
   
         Eigen::SelfAdjointEigenSolver<Mat> ces;
         ces.compute(K);
   
         wK_ = ces.eigenvalues();  
         uK_ = ces.eigenvectors(); 
         uKdag_ = ces.eigenvectors().adjoint(); 

         //std::cout << "K:\n" << K << std::endl; 
         //std::cout << "wK_:\n" << wK_ << std::endl; 
         //std::cout << "uK_:\n" << uK_ << std::endl; 
 
         //std::cout << "U*Udag:\n" << uK_ * uKdag_ << std::endl; 
        
         init_without_vertex(); 
        }

        void init_without_vertex(){

           U_ = Mat::Identity(ns_, ns_); 
           D_ = Mat::Zero(ns_, ns_); 
           V_ = Mat::Identity(ns_, ns_);

           //initially gtau_ is noninteracting gf = 1./(1+ exp(-E*beta))
           for(site_type l=0; l<ns_; ++l) 
               D_(l, l) = wK_(l)>0. ? 1./(1.+exp(-beta_*wK_(l))) : exp(beta_*wK_(l))/(1.+exp(beta_*wK_(l))) ; // it is samething, to avoid overflow 
           
           //std::cout << "initially:" << std::endl; 
           //std::cout << "D_:\n" << D_<< std::endl; 
           
           //fill in storage
           //since initially there is no vertices U and V are all diagonal 
           for (unsigned ib=0; ib< Rstorage_.size(); ++ib) {
               Vec v(ns_); 
               for (site_type l=0; l<ns_; ++l)
                   v(l) = exp(-(itime2time(ib*blocksize_)) *wK_(l)); 
               Rstorage_[ib] = boost::make_tuple(Mat::Identity(ns_, ns_), v, Mat::Identity(ns_, ns_)); 
           }
           
           for (unsigned ib=0; ib< Lstorage_.size(); ++ib) {
               Vec v(ns_); 
               for (site_type l=0; l<ns_; ++l)
                   v(l) = exp(-(beta_ -itime2time((ib+1)*blocksize_)) *wK_(l)); 
               Lstorage_[ib] = boost::make_tuple(Mat::Identity(ns_, ns_), v, Mat::Identity(ns_, ns_)); 
           }
        }

        void rebuild(const tlist_type& tlist, vlist_type& vlist){
            
            Mat U, D, V; 
            boost::tie(U, D, V) = Gstable(itau_, tlist, vlist); // from scratch 

            double max_diff = ((U*D*V - U_*D_*V_).cwiseAbs()).maxCoeff(); 
            if(max_diff > 1.e-6){
              std::cout<<"WARNING: roundoff errors " <<max_diff << std::endl;

              //std::cout << "in rebuild:" << std::endl; 
              //std::cout << "U_:\n" << U_ << std::endl; 
              //std::cout << "D_:\n" << D_ << std::endl; 
              //std::cout << "V_:\n" << V_ << std::endl; 

              //std::cout << "U:\n" << U << std::endl; 
              //std::cout << "D:\n" << D << std::endl; 
              //std::cout << "V:\n" << V << std::endl; 
          
              //std::cout << "U_*D_*V_:\n" << U_*D_*V_ << std::endl; 
              //std::cout << "U*D*V:\n" << U*D*V << std::endl; 
              //std::cout << "diff:\n" <<U_*D_*V_-U*D*V << std::endl; 

              //std::cout << "tlist: "; 
              //std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
              //std::cout << std::endl; 
            }
           
            U_ = U;  
            D_ = D; 
            V_ = V;
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
        //but does not affect the storage 
        void update(const site_type si, const site_type sj, const double gij){

             //std::cout << "si, sj, gij: "  << si << " " << sj << " " << gij << std::endl; 
             //std::cout << "uKdag_.col(si):\n" << uKdag_.col(si) << std::endl; 
             //std::cout << "uKdag_.col(sj):\n" << uKdag_.col(sj) << std::endl; 
             //std::cout << "uK_.row(si):\n" << uK_.row(si) << std::endl; 
             //std::cout << "uK_.row(sj):\n" << uK_.row(sj) << std::endl; 

             //update D_
             D_ -=  ( (D_* (V_* uKdag_.col(sj))) * ((uK_.row(si)*U_)*D_)
                    + (D_* (V_* uKdag_.col(si))) * ((uK_.row(sj)*U_)*D_)
                    )/gij; 
             
             // * U^\dagger V U from right 
             V_ -= 2.* (V_* uKdag_.col(si)) * uK_.row(si) + 2.* (V_*uKdag_.col(sj)) * uK_.row(sj); // this thing is fixed 

             if (itau_%blocksize_==0){ //special treatment when update the block starting point 
                  //std::cout << "update: special treatment because update on the block boundary" << std::endl; 
                  Vprop(si, sj, "L",  Rstorage_[itau_/blocksize_].get<0>());// update U in Rstorage 
             }

             if (update_refresh_counter_ < update_refresh_period_){
                    ++update_refresh_counter_; 
             }else{

                    Eigen::JacobiSVD<Mat> svd(D_, Eigen::ComputeThinU | Eigen::ComputeThinV); 
                    U_ = U_*svd.matrixU(); 
                    D_ = svd.singularValues().asDiagonal(); 
                    V_ = svd.matrixV().adjoint()*V_;
                    
                    update_refresh_counter_ =0; 
             }

             //std::cout << "in update:" << std::endl; 
             //std::cout << "U_:\n" << U_ << std::endl; 
             //std::cout << "D_:\n" << D_ << std::endl; 
             //std::cout << "V_:\n" << V_ << std::endl; 

        }
        
        //wrap does not change D_  
        void wrap(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {
            if (itau == itau_) 
                return; 
            
             itime_type b = itau/blocksize_; //new block index 
             itime_type b_ = itau_/blocksize_; //old block index 

             //std::cout << "b, b_, itau, itau_, blocksize: "  << b << " " << b_ << " " << itau << " " << itau_ << " " << blocksize_ << std::endl; 

             bool refresh = false; 
             if (wrap_refresh_counter_ < wrap_refresh_period_){
                 ++wrap_refresh_counter_; 
             }else{

                 refresh = true; 
                 wrap_refresh_counter_ = 0; 
             }

            //wrap Green's function 
            if (itau >= itau_) {
                // B G B^{-1}
                propagator1(-1, itau, itau_, tlist, vlist, U_);  // B(tau1) ... B(tau2) *U_  

                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(U_*D_, Eigen::ComputeThinU | Eigen::ComputeThinV); 
                  U_ = svd.matrixU();
                  D_ = svd.singularValues().asDiagonal();
                  V_ = svd.matrixV().adjoint()*V_;
                }

                propagator1(1, itau, itau_, tlist, vlist, V_); // V_ * B^{-1}(tau2) ... B^{-1}(tau1)
                
                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(D_*V_, Eigen::ComputeThinU | Eigen::ComputeThinV); 
                  U_ = U_*svd.matrixU();
                  D_ = svd.singularValues().asDiagonal();
                  V_ = svd.matrixV().adjoint();
                }

                itau_ = itau; 

            }else{

                // B^{-1} G B 
                propagator2(1, itau_, itau, tlist, vlist, U_); //  B^{-1}(tau2) ... B^{-1}(tau1) * U_
                    
                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(U_*D_, Eigen::ComputeThinU | Eigen::ComputeThinV); 
                  U_ = svd.matrixU();
                  D_ = svd.singularValues().asDiagonal();
                  V_ = svd.matrixV().adjoint()*V_;
                }

                propagator2(-1, itau_, itau, tlist, vlist, V_);   //  V_ * B(tau1) ... B(tau2)
                    
                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(D_*V_, Eigen::ComputeThinU | Eigen::ComputeThinV); 
                  U_ = U_*svd.matrixU();
                  D_ = svd.singularValues().asDiagonal();
                  V_ = svd.matrixV().adjoint();
                }

                itau_ = itau; 
            }

            //when we wrap to a new block we need to update storage 
            if (b > b_){// move to a larger block on the left  
                assert(b-b_ ==1) ; 
                    
                Mat U,V;  
                Vec D; 
                boost::tie(U, D, V) = Rstorage_[b_];
                propagator1(-1, b*blocksize_, b_*blocksize_, tlist, vlist, U);

                if (refresh){
                   Eigen::JacobiSVD<Mat> svd(U*D.asDiagonal(), Eigen::ComputeThinU | Eigen::ComputeThinV); 
                   U = svd.matrixU();
                   D = svd.singularValues();
                   V = svd.matrixV().adjoint()*V;
                }

                Rstorage_[b] = boost::make_tuple(U, D, V); 

            }else if (b< b_){// move to smaller block 
                assert(b_-b ==1); 
                
                Mat U, V;  
                Vec D; 
                boost::tie(U, D, V) = Lstorage_[b_];
                propagator2(-1, (b_+1)*blocksize_, b_*blocksize_, tlist, vlist, V);

                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(D.asDiagonal()*V, Eigen::ComputeThinU | Eigen::ComputeThinV); 
                  U = U*svd.matrixU();
                  D = svd.singularValues();
                  V = svd.matrixV().adjoint();
                }

                Lstorage_[b] = boost::make_tuple(U, D, V); 
            }
        }

         /*
         //equal time Green's function at tau 
         //for test only 
        boost::tuple<Mat, Mat, Mat> G(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {// this is very expansive because of inverse
           //Mat res = Mat::Identity(ns_, ns_) + B(itau, 0, tlist, vlist) * B(itime_max, itau, tlist, vlist); 
           //Mat res = Mat::Identity(ns_, ns_) + B_tau_0 * B_beta_tau; 
       
            Mat B_tau_0 = Mat::Identity(ns_, ns_); 
            propagator1(-1, itau, 0, tlist, vlist, B_tau_0);
            
            Mat B_beta_tau =  Mat::Identity(ns_, ns_);
            propagator2(-1, itime_max, itau, tlist, vlist, B_beta_tau);
            
            Mat res = Mat::Identity(ns_, ns_) + B_tau_0 * B_beta_tau; 

            Eigen::JacobiSVD<Mat> svd(res.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV); 

            return boost::tie(svd.matrixU(), svd.singularValues().asDiagonal(), svd.matrixV().adjoint() ) ; 
         }
        */
        
         //we return U,D,V but not G 
         boost::tuple<Mat, Mat, Mat> Gstable(const itime_type itau, const tlist_type& tlist, vlist_type& vlist)  const {

           itime_type b = itau/blocksize_; //block index 
           //std::cout << "itau, block: " << itau << " " << b << std::endl; 

           //B_tau_0 = U1*D1*V1
           Mat U1, V1;  
           Vec D1; 
           boost::tie(U1, D1, V1) = Rstorage_[b]; 
           propagator1(-1, itau, b*blocksize_, tlist, vlist, U1);
        
           //B_beta_tau = U2*D2*V2
           Mat U2, V2;  
           Vec D2;
           boost::tie(U2, D2, V2) = Lstorage_[b]; 
           propagator2(-1, (b+1)*blocksize_, itau, tlist, vlist, V2);

           //std::cout << "in Gstable:" << std::endl; 

           //std::cout << "U1:\n" << U1 << std::endl; 
           //std::cout << "D1:\n" << D1 << std::endl; 
           //std::cout << "V1:\n" << V1 << std::endl; 

           //std::cout << "U2:\n" << U2 << std::endl; 
           //std::cout << "D2:\n" << D2 << std::endl; 
           //std::cout << "V2:\n" << V2 << std::endl; 

        
           /*(3)
           Eigen::JacobiSVD<Mat> svd((D1.asDiagonal()*V1)*(U2*D2.asDiagonal()), Eigen::ComputeThinU | Eigen::ComputeThinV); 
           Mat U, D, V;  
           U = U1 * svd.matrixU(); 
           D = svd.singularValues().asDiagonal(); 
           V = svd.matrixV().adjoint()*V2;

           Mat res= U.inverse()*V.inverse() + D;

           Eigen::JacobiSVD<Mat> svd2(res, Eigen::ComputeThinU | Eigen::ComputeThinV); 
           D =  svd2.singularValues().asDiagonal();

           return boost::make_tuple((svd2.matrixV().adjoint()*V).inverse() ,  D.inverse() ,  (U*svd2.matrixU()).inverse()); 
           */

           /*(1)
           Mat res = Mat::Identity(ns_, ns_) + (U1*D1.asDiagonal()*V1) * (U2*D2.asDiagonal()*V2); 
           Eigen::JacobiSVD<Mat> svd(res.inverse(), Eigen::ComputeThinU | Eigen::ComputeThinV); 
           return boost::tie(svd.matrixU(), svd.singularValues().asDiagonal(), svd.matrixV().adjoint()); 
           */
          
           //(2)
           Mat res= (V2*U1).inverse() + D1.asDiagonal() * (V1* U2) * D2.asDiagonal();
           Eigen::JacobiSVD<Mat> svd(res, Eigen::ComputeThinU | Eigen::ComputeThinV); 
           
           Mat U, D, V;   
           U = svd.matrixU(); 
           D = svd.singularValues().asDiagonal(); 
           V = svd.matrixV().adjoint();

           //std::cout << "in Gstable:" << std::endl; 
           //std::cout << "U:\n" << U << std::endl; 
           //std::cout << "D:\n" << D << std::endl; 
           //std::cout << "V:\n" << V << std::endl; 

           return boost::make_tuple((V*V2).inverse() ,  D.inverse() ,  (U1*U).inverse()); 

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
         // Btau(2) does not contain vertex contribution  
         void propagator1(const int sign, const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist, Mat& A)const  { 
     
             assert(itau1>=itau2); 

             std::string side = (sign ==-1) ? "L" : "R"; 
            
             //since we always point to the right of a vertex, we push downwards the lower boundary 
             //this is different with block convention (,]  
             tlist_type::const_iterator lower, upper; 
             lower = std::lower_bound (tlist.begin(), tlist.end(), itau2, std::less_equal<itime_type>()); //equal is exclude
             upper = std::upper_bound (tlist.begin(), tlist.end(), itau1); 

             //std::cout << "props1: itau1, itau2, vertices at " << itau1 << " " << itau2 << std::endl; 
             //std::copy(lower, upper, std::ostream_iterator<itime_type>(std::cout, " "));
             //std::cout << std::endl;  

             //upper > tau1 > lower > tau2 
             if (lower == upper ) {// there is no vertex in between tau1 and tau2 
                 Kprop(sign, itau1 - itau2, side, A); 

             }else{

                 Kprop(sign, *lower - itau2, side, A);
                 //std::cout << "prop1:Kprop " <<   *lower  << " " << itau2  << " " << *lower - itau2 << std::endl; 
                 for (tlist_type::const_iterator it1 =lower, it2 =++lower; it1!=upper; ++it1, ++it2) {
                    
                     itime_type itau = *it1; 
                     Vprop(vlist[itau][0], vlist[itau][1], side, A); 
                     //std::cout << "prop1:act vertex at " << itau  << std::endl; 
                 
                     itime_type ditau = (it2 ==upper) ? itau1 - itau: *it2 - itau; 
                     Kprop(sign, ditau , side, A); 
                     //std::cout << "prop1:Kprop " <<   ((it2 ==upper) ? itau1: *it2)  << " " << itau  << " " << ditau << std::endl; 
                 }
                 //std::cout << "##############" << std::endl; 
             }
         }

         // it can do A* B(tau1)... B(tau2)  when sign = -1
         // or        B(tau2)^{-1} ... B(tau1)^{-1}* A when sign = 1 
         void propagator2(const int sign, const itime_type itau1, const itime_type itau2, const tlist_type& tlist, vlist_type& vlist, Mat& A) const { 
     
             assert(itau1>=itau2); 

             std::string side = (sign ==1) ? "L" : "R"; 
     
             tlist_type::const_iterator lower, upper; 
             lower = std::lower_bound (tlist.begin(), tlist.end(), itau2, std::less_equal<itime_type>()); //equal is exclude
             upper = std::upper_bound (tlist.begin(), tlist.end(), itau1);
                
             //std::cout << "props2: itau1, itau2, vertices at " << itau1 << " " << itau2 << std::endl; 
             //std::copy(lower, upper, std::ostream_iterator<itime_type>(std::cout, " "));
             //std::cout << std::endl;  

             //upper > tau1 > lower > tau2 
             if (lower == upper ) {// there is no vertex in between tau1 and tau2 
                 Kprop(sign, itau1 - itau2, side, A); 

             }else{

                 Kprop(sign, itau1 - *(--upper), side, A);
                 //std::cout << "prop2:Kprop " <<  itau1 << " " << *upper  << " " << itau1 - *upper << std::endl; 

                 for (tlist_type::const_iterator it1 =upper, it2 =--upper; it1!=lower; --it1, --it2) {
                    
                     itime_type itau = *it1; 
                     Vprop(vlist[itau][0], vlist[itau][1], side, A); 
                     //std::cout << "prop2:act vertex at " << itau  << std::endl; 
                 
                     itime_type ditau = itau - *it2; 
                     Kprop(sign, ditau , side, A); 
                     //std::cout << "prop2:Kprop " <<  itau << " " << *it2  << " " << ditau << std::endl; 
                 }
                 {
                     //the last step by hand (it1 = lower, it2 point to tau2) 
                     itime_type itau = *lower; 
                     Vprop(vlist[itau][0], vlist[itau][1], side, A); 
                     //std::cout << "prop2:act vertex at " << itau  << std::endl; 
             
                     itime_type ditau = itau - itau2; 
                     Kprop(sign, ditau , side, A); 
                     //std::cout << "prop2:Kprop " <<  itau << " " <<  itau2 << " " << ditau << std::endl; 
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
            return timestep_*boost::lexical_cast<double>(itau);  
        }

    private:
        const site_type ns_; 
        const time_type beta_; 
        const time_type timestep_; 

        //eigen value and vectors of K 
        Vec wK_; 
        Mat uK_; 
        Mat uKdag_; 

        itime_type itau_; 

        //gtau = U_*D_*V_
        Mat U_, D_, V_; 
        
        //blocksize is used in wrap,when determine the block index 
        const itime_type blocksize_;

        //counter and period for refreshing 
        unsigned update_refresh_counter_; 
        unsigned update_refresh_period_; 

        unsigned wrap_refresh_counter_; 
        unsigned wrap_refresh_period_; 

        //storage
        std::vector<boost::tuple<Mat, Vec, Mat> > Lstorage_;
        std::vector<boost::tuple<Mat, Vec, Mat> > Rstorage_;

};

#endif
