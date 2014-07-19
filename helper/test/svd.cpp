//g++ svd.cpp -I/users/lewang/Libs/include/eigen3/
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using Eigen::JacobiSVD;

using std::cout;
using std::endl;


int main(){

    MatrixXd m = MatrixXd::Random(500,500);

    //cout << "Here is the matrix m:" << endl << m << endl;
    JacobiSVD<MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

    return 0; 

    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
   
    cout << "U*D*V^T"  << endl; 
    cout << svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().adjoint() << endl;  
    cout << svd.matrixV() *  svd.matrixV().adjoint() << endl;  
    cout << svd.matrixU() *  svd.matrixU().adjoint() << endl;  

    return 0; 
}
