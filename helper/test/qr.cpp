#include <Eigen/Core>
#include <Eigen/QR>
#include <iostream>

using Eigen::MatrixXf;
using Eigen::HouseholderQR;

int main(){

    MatrixXf A(MatrixXf::Random(4,4)), Q, R;

    A.setRandom();
    HouseholderQR<MatrixXf> qr(A);
    Q = qr.householderQ(); 
    R = qr.matrixQR().triangularView<Eigen::Upper>();

    std::cout << "Q:\n" << Q << std::endl; 
    std::cout << "R:\n" << R << std::endl; 

    std::cout << "A:\n" << A << std::endl; 
    std::cout << "Q*R:\n" << Q*R << std::endl; 

}

