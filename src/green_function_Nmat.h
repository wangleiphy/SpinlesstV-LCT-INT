#ifndef GREEN_FUNCTION_NMAT_H
#define GREEN_FUNCTION_NMAT_H

#include "types.h"
#include <boost/tuple/tuple.hpp>
#include <iterator>
#include <iostream>

class Green_function{

    public:

        Green_function(const Mat& K, const time_type timestep, const unsigned nblock, const itime_type blocksize, const unsigned wrap_refresh_period)
        :ns_(K.rows())
        ,np_(ns_/2)// half filled 
        ,timestep_(timestep)
        ,itau_(0)
        ,R_(Mat::Identity(ns_, np_))// these are the regularied version 
        ,N_(Mat::Identity(np_, np_))
        ,L_(Mat::Identity(np_, ns_))
        ,blocksize_(blocksize)
        ,wrap_refresh_counter1_(0)
        ,wrap_refresh_counter2_(0)
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

          //fill in storage
           //since initially there is no vertices U and V are all diagonal 
           for (unsigned ib=0; ib< Rstorage_.size(); ++ib) {
               Rstorage_[ib] = Mat::Identity(ns_, np_); // U^{dagger} *P and assume P is the eigen state of K 
           }
           
           for (unsigned ib=0; ib< Lstorage_.size(); ++ib) {
               Lstorage_[ib] = Mat::Identity(np_, ns_); 
           }
        }

        void rebuild(const tlist_type& tlist, vlist_type& vlist){

            Mat R, N, L; 
            boost::tie(R, N, L) = stablization(itau_, tlist, vlist); // from scratch 

            double max_diff = ((R*N*L - R_*N_*L_).cwiseAbs()).maxCoeff(); 
            if(max_diff > 1.e-6){
              std::cout<<"WARNING: roundoff errors " <<max_diff << std::endl;

              //std::cout << "in rebuild:" << std::endl; 
         
              //std::cout << "gtau:\n" << R*N*L << std::endl; 
              //std::cout << "gtau_:\n" << R_*N_*L_ << std::endl; 
              //std::cout << "diff:\n" << R*N*L - R_*N_*L_ << std::endl; 

              //std::cout << "tlist: "; 
              //std::copy(tlist.begin(), tlist.end(), std::ostream_iterator<itime_type>(std::cout, " "));
              //std::cout << std::endl; 
              //abort(); 
            }

            R_ = R;
            N_ = N;
            L_ = L; 
        }
        
        Mat gtau() const {
            return Mat::Identity(ns_, ns_) - R_*N_*L_; 
        }

        time_type tau() const {
            return itime2time(itau_); 
        }

        double gij(const site_type si, const site_type sj)const {
            double delta = (si==sj)? 1.0 : 0.0; 
            return delta -(uK_.row(si) *R_)*  N_ * (L_*uKdag_.col(sj));  
        }

        Vec denmat(const site_type si)const {
            Vec res = Vec::Zero(ns_); 
            res(si) = 1.0; 
            return res  -  uK_ *( R_* (N_*  (L_ * uKdag_.col(si))));
        }
            
        //update changes D_ and V_ 
        //but does not affect the storage 
        void update(const site_type si, const site_type sj, const double gij){

             //std::cout << "si, sj, gij: "  << si << " " << sj << " " << gij << std::endl; 
             //std::cout << "uKdag_.col(si):\n" << uKdag_.col(si) << std::endl; 
             //std::cout << "uKdag_.col(sj):\n" << uKdag_.col(sj) << std::endl; 
             //std::cout << "uK_.row(si):\n" << uK_.row(si) << std::endl; 
             //std::cout << "uK_.row(sj):\n" << uK_.row(sj) << std::endl; 

             //std::cout << "Lshape: " << L_.rows() << " " << L_.cols() << std::endl;  
             //std::cout << "Nshape: " << N_.rows() << " " << N_.cols() << std::endl;  
             //std::cout << "Rshape: " << R_.rows() << " " << R_.cols() << std::endl;  

             //update N_
             N_ +=  ((N_* (L_* uKdag_.col(sj))) * ((uK_.row(si)*R_)*N_)
                     +(N_* (L_* uKdag_.col(si))) * ((uK_.row(sj)*R_)*N_)
                    )/gij; 

             //update R_
             Vprop(si, sj, "L",  R_);

             if (itau_%blocksize_==0){ //special treatment when update the block starting point 
                  //std::cout << "update: special treatment because update on the block boundary" << std::endl; 
                  Rstorage_[itau_/blocksize_] = R_; 
             }

             //std::cout << "in update:" << std::endl; 
             //std::cout << "U_:\n" << U_ << std::endl; 
             //std::cout << "D_:\n" << D_ << std::endl; 
             //std::cout << "V_:\n" << V_ << std::endl; 
        }
        
        //wrap does not change N_ 
        void wrap(const itime_type itau, const tlist_type& tlist, vlist_type& vlist) {
            if (itau == itau_) 
                return; 

            
             itime_type b = itau/blocksize_; //new block index 
             itime_type b_ = itau_/blocksize_; //old block index 

             //std::cout << "b, b_, itau, itau_, blocksize: "  << b << " " << b_ << " " << itau << " " << itau_ << " " << blocksize_ << std::endl; 

            //refresh every wrap_refresh_period_ steps 
             bool refresh = false; 
             if (wrap_refresh_counter1_ < wrap_refresh_period_){
                 ++wrap_refresh_counter1_; 
             }else{
                 refresh = true; 
                 wrap_refresh_counter1_ = 0; 
             }

            //wrap Green's function 
            if (itau >= itau_) {
                // B G B^{-1}
                propagator1(-1, itau, itau_, tlist, vlist, R_);  // B(tau1) ... B(tau2) *R_  
                
                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(R_, Eigen::ComputeThinU); 
                  R_ = svd.matrixU();
                }

                propagator1(1, itau, itau_, tlist, vlist, L_); // L_ * B^{-1}(tau2) ... B^{-1}(tau1)

                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(L_, Eigen::ComputeThinV); 
                  L_ = svd.matrixV().adjoint();
                  N_ = (L_*R_).inverse();   // since we have changed L_ and R_, we need to change N_ as well 
                }

                itau_ = itau; 

            }else{

                // B^{-1} G B 
                propagator2(1, itau_, itau, tlist, vlist, R_); //  B^{-1}(tau2) ... B^{-1}(tau1) * R_
                   
                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(R_, Eigen::ComputeThinU); 
                  R_ = svd.matrixU();
                }

                propagator2(-1, itau_, itau, tlist, vlist, L_);   //  L_ * B(tau1) ... B(tau2)

                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(L_, Eigen::ComputeThinV); 
                  L_ = svd.matrixV().adjoint();
                  N_ = (L_*R_).inverse();
                }

                itau_ = itau; 
            }
            
            //refresh every wrap_refresh_period_ blocks  
             refresh = false; 
             if (wrap_refresh_counter2_ < wrap_refresh_period_){
                 ++wrap_refresh_counter2_; 
             }else{
                 refresh = true; 
                 wrap_refresh_counter2_ = 0; 
             }

            //when we wrap to a new block we need to update storage 
            if (b > b_){// move to a larger block on the left  
                assert(b-b_ ==1) ; 
                    
                Mat UR= Rstorage_[b_];
                propagator1(-1, b*blocksize_, b_*blocksize_, tlist, vlist, UR);

                if (refresh){
                   Eigen::JacobiSVD<Mat> svd(UR, Eigen::ComputeThinU); 
                   UR = svd.matrixU();
                }

                Rstorage_[b] = UR; 

            }else if (b< b_){// move to smaller block 
                assert(b_-b ==1); 
                
                Mat VL = Lstorage_[b_];
                propagator2(-1, (b_+1)*blocksize_, b_*blocksize_, tlist, vlist, VL);

                if (refresh){
                  Eigen::JacobiSVD<Mat> svd(VL, Eigen::ComputeThinV); 
                  VL = svd.matrixV().adjoint();
                }

                Lstorage_[b] = VL; 
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
        
         
        boost::tuple<Mat, Mat, Mat> stablization(const itime_type itau, const tlist_type& tlist, vlist_type& vlist)  const {

           itime_type b = itau/blocksize_; //block index 
           //std::cout << "itau, block: " << itau << " " << b << std::endl; 

           Mat R = Rstorage_[b]; 
           propagator1(-1, itau, b*blocksize_, tlist, vlist, R);
        
           Mat L = Lstorage_[b]; 
           propagator2(-1, (b+1)*blocksize_, itau, tlist, vlist, L);

           return boost::make_tuple(R, (L*R).inverse(), L); 
         }
                
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
        const site_type np_; 
        const time_type timestep_; 

        //eigen value and vectors of K 
        Vec wK_; 
        Mat uK_; 
        Mat uKdag_; 

        itime_type itau_; 
        
        //gtau = 1- R_ * N_ * L_; 
        //N_ = (L_*R_)^{-1}
        Mat R_, N_, L_; 
        
        //blocksize is used in wrap,when determine the block index 
        const itime_type blocksize_;

        //counter and period for refreshing 
        //unsigned update_refresh_counter_; 
        //unsigned update_refresh_period_; 

        unsigned wrap_refresh_counter1_; 
        unsigned wrap_refresh_counter2_; 

        unsigned wrap_refresh_period_; 

        //storage
        std::vector<Mat> Lstorage_;
        std::vector<Mat> Rstorage_;

};

#endif