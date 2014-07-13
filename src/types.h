#ifndef TYPES_H
#define TYPES_H
#include <vector>
#include <map>
#include <set>
#include <Eigen/Dense>

typedef unsigned site_type;  
typedef unsigned itime_type;  

typedef double time_type; 
typedef std::set<itime_type> tlist_type; 
typedef std::map<itime_type, std::vector<site_type> > vlist_type;

typedef Eigen::MatrixXd Mat;  

static const itime_type itime_max = std::numeric_limits<itime_type>::max(); 

#endif /*TYPES_H_*/
