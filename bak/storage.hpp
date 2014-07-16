#ifndef STORAGE_H
#define STORAGE_H

#include "types.h"

class Storage{
   
    public: 
        Storage(unsigned nblock)
        :U_(nblock) 
        ,D_(nblock) 
        ,V_(nblock) 
        {}

        const Mat& U(unsigned i)const { return U_[i];}
        const Mat& D(unsigned i)const { return D_[i];}
        const Mat& V(unsigned i)const { return V_[i];}

        Mat& U(unsigned i){ return U_[i];}
        Mat& D(unsigned i){ return D_[i];}
        Mat& V(unsigned i){ return V_[i];}
       
    private:
        std::vector<Mat> U_; 
        std::vector<Mat> D_; 
        std::vector<Mat> V_; 
}; 

#endif 
