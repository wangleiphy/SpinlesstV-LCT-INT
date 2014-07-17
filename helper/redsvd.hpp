/* 
 *  Copyright (c) 2010 Daisuke Okanohara
 * 
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 * 
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#ifndef REDSVD_HPP__
#define REDSVD_HPP__

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "util.hpp"

namespace REDSVD {

class RedSVD {
public:
  RedSVD(){}

  template <class Mat>
  RedSVD(Mat& A){
    int r = (A.rows() < A.cols()) ? A.rows() : A.cols();
    run(A, r);
  }

  template <class Mat>
  RedSVD(Mat& A, const int rank){
    run(A, rank);
  }

  template <class Mat>
  void run(Mat& A, const int rank){
    if (A.cols() == 0 || A.rows() == 0) return;
    int r = (rank < A.cols()) ? rank : A.cols();
    r = (r < A.rows()) ? r : A.rows();
    
    // Gaussian Random Matrix for A^T
    Eigen::MatrixXd O(A.rows(), r);
    Util::sampleGaussianMat(O);
    
    // Compute Sample Matrix of A^T
    Eigen::MatrixXd Y = A.transpose() * O;
    
    // Orthonormalize Y
    Util::processGramSchmidt(Y);

    // Range(B) = Range(A^T)
    Eigen::MatrixXd B = A * Y;
    
    // Gaussian Random Matrix
    Eigen::MatrixXd P(B.cols(), r);
    Util::sampleGaussianMat(P);
    
    // Compute Sample Matrix of B
    Eigen::MatrixXd Z = B * P;
    
    // Orthonormalize Z
    Util::processGramSchmidt(Z);
    
    // Range(C) = Range(B)
    Eigen::MatrixXd C = Z.transpose() * B; 
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svdOfC(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    // C = USV^T
    // A = Z * U * S * V^T * Y^T()
    matU_ = Z * svdOfC.matrixU();
    matS_ = svdOfC.singularValues();
    matV_ = Y * svdOfC.matrixV();
  }
  
  const Eigen::MatrixXd& matrixU() const {
    return matU_;
  }

  const Eigen::VectorXd& singularValues() const {
    return matS_;
  }

  const Eigen::MatrixXd& matrixV() const {
    return matV_;
  }

private:
  Eigen::MatrixXd matU_;
  Eigen::VectorXd matS_;
  Eigen::MatrixXd matV_;
};

class RedSymEigen {
public:
  RedSymEigen(){}

  template <class Mat>
  RedSymEigen(Mat& A, const int rank){
    run(A, rank);
  }  

  template <class Mat>
  void run(Mat& A, const int rank){
    if (A.cols() == 0 || A.rows() == 0) return;
    int r = (rank < A.cols()) ? rank : A.cols();
    r = (r < A.rows()) ? r : A.rows();
    
    // Gaussian Random Matrix
    Eigen::MatrixXd O(A.rows(), r);
    Util::sampleGaussianMat(O);
    
    // Compute Sample Matrix of A
    Eigen::MatrixXd Y = A.transpose() * O;
    
    // Orthonormalize Y
    Util::processGramSchmidt(Y);

    Eigen::MatrixXd B = Y.transpose() * A * Y;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenOfB(B);
    
    eigenValues_ = eigenOfB.eigenvalues();
    eigenVectors_ = Y * eigenOfB.eigenvectors();
  }
  
  const Eigen::MatrixXd& eigenVectors() const {
    return eigenVectors_;
  }

  const Eigen::VectorXd& eigenValues() const {
    return eigenValues_;
  }

private:
  Eigen::VectorXd eigenValues_;
  Eigen::MatrixXd eigenVectors_;
};

class RedPCA {
public:
  RedPCA(){}

  template <class Mat>
  RedPCA(const Mat& A, const int rank) {
    run(A, rank);
  }

  template <class Mat> 
  void run(const Mat& A, const int rank) {
    RedSVD redsvd;
    redsvd.run(A, rank);
    const Eigen::VectorXd& S = redsvd.singularValues();
    principalComponents_ = redsvd.matrixV();
    scores_              = redsvd.matrixU() * S.asDiagonal();
  }

  const Eigen::MatrixXd& principalComponents() const {
    return principalComponents_;
  }

  const Eigen::MatrixXd& scores() const {
    return scores_;
  }

 private:
  Eigen::MatrixXd principalComponents_;
  Eigen::MatrixXd scores_;
};

}

#endif // REDSVD_HPP__
