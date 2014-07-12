#ifndef MMATRIX_HPP
#define MMATRIX_HPP

#include "operator.hpp"
#include <Eigen/Dense>

class m_matrix
{
public:
    typedef Eigen::MatrixXd matrix_t ; 

    m_matrix()
        :num_vertices_(0)
        ,creators_()
        ,matrix_()
    {}
    
    void clear(){
      num_vertices_ = 0; 
      creators_.clear();
      matrix_.resize(0,0);  
    }

    matrix_t& matrix() { return matrix_;}
    matrix_t const &matrix() const { return matrix_;}

    std::vector<creator> &creators(){ return creators_;}
    const std::vector<creator> &creators() const{ return creators_;}

    unsigned int & num_vertices(){ return num_vertices_;}
    const unsigned int & num_vertices()const{ return num_vertices_;}

private:
    unsigned int num_vertices_; 

    std::vector<creator> creators_;         //an array of creation operators c_dagger corresponding to the row of the matrix

    matrix_t matrix_;

};

#endif 
