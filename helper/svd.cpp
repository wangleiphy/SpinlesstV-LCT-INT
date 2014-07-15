#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>

using Eigen::MatrixXf;
using Eigen::VectorXf;

using Eigen::JacobiSVD;

using std::cout;
using std::endl;


int main(){

    MatrixXf m = MatrixXf::Random(4,4);
    cout << "Here is the matrix m:" << endl << m << endl;
    JacobiSVD<MatrixXf> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
   
    cout << "U*D*V^T"  << endl; 
    cout << svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().adjoint() << endl;  

    return 0; 
}
