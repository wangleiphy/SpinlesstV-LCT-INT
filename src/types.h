#ifndef TYPES_H
#define TYPES_H
#include <vector>
#include <map>
#include <set>
//#include <limits>
#include <Eigen/Dense>

typedef unsigned site_type;  
typedef unsigned itime_type;  

typedef double time_type; 
typedef std::set<itime_type> tlist_type; 
typedef std::map<itime_type, std::vector<site_type> > vlist_type;

typedef Eigen::MatrixXd Mat;  
typedef Eigen::VectorXd Vec; 

typedef std::map<unsigned, std::vector<site_type> > DistanceMap;

#endif /*TYPES_H_*/
