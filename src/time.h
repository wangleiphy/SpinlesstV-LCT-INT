#ifndef TIME_H
#define TIME_H

#include "types.h"

time_type index2time(const index_type itime){
    return beta * itime/std::numeric_limits<unsigned long>::max(); 
}

#endif 
