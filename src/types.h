#ifndef TYPES_H
#define TYPES_H
#include <vector>
#include <map>
#include <set>
#include <Eigen/Dense>

typedef unsigned site_type;  
typedef double time_type;  

typedef std::set<time_type> tlist_type; 
typedef std::map<time_type, std::vector<site_type> > vlist_type;

typedef Eigen::MatrixXd Mat;  

#endif /*TYPES_H_*/
