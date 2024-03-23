#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <math.h>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
Eigen::MatrixXd K1_cov_est(NumericVector x, NumericVector dim, const int k) {
  // NumericVector x = clone(xi); (for debug, no need now)
  // k is mode starting from 2 (using the way that R counts)
  // since time mode is mode 1;
  // so for K=2, we have k=2 ;

  int dim1 = dim[0];
  int dim2 = dim[1];
  int dimkm1 = dim[k-1];
  
  Eigen::Tensor<double, 2> tensor(dim1, dim2);
  
  for (int i=0; i < dim1; i++){
    for (int j=0; j < dim2; j++){
      tensor(i,j) = x(i + j*dim1);
    }
  }
  // finish loading the tensor by now
  
  
  
  Eigen::MatrixXd weight_cov = Eigen::MatrixXd::Constant(dimkm1, dimkm1, 0.0);
  
  int rows = dimkm1;
  long cols = tensor.size() / rows;
  cols /= dim1;
  
  
  // create unfolding first, for each t
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> matrixList(dim1);
  
  for (int t = 0; t < dim1; t++) {
    // Create an Eigen Matrix to store the unfolded data
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> k_mat(rows, cols);
    
    // Unfold the tensor into a matrix
    for (int u1 = 0; u1 < rows; u1++) {
      for (long u2 = 0; u2 < cols; u2++) {
        long index = u1 * cols + u2;
        Eigen::array<Eigen::Index, 2> indices;
        indices[0] = t;
        indices[k-1] = u1;
        for (int w = 1; w < 2; w++) {
          if ( w != (k-1)) {
            indices[w] = int( fmod(index, dim[w]) );
            index /= dim[w];
          }
        }
        
        k_mat(u1, u2) = tensor(indices);
      }
    }
    
    // matrixList.push_back(k_mat);
    matrixList[t] = k_mat;
  }
  
  
  // create Boolean unfolding to check NA first, for each t
  std::vector<Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>> observeList(dim1);
  for (int t = 0; t < dim1; t++) {
    // Create an Eigen Matrix to store the unfolded data
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> observe_mat(rows, cols);
    
    for (int i = 0; i < rows; i++){
      for (long j = 0; j < cols; j++){
        observe_mat(i,j) = NumericVector::is_na( matrixList[t](i,j) );
      }
    }
    
    
    // observeList.push_back( observe_mat );
    observeList[t] = observe_mat;
  }
  
  
  
  
  for (int i = 0; i < dimkm1; i++) {
    //for (int j = i + 1; j < weight_cov.cols(); j++) {
    for (int j = i; j < dimkm1; j++) {
      long double cov_acc = 0;
      
      for (long c = 0; c < cols; c++){
        int obs_t = 0;
        long double c_cov = 0;
        
        for (int t = 0; t < dim1; t++){
          
          
          if ( !( observeList[t](i,c) ) && !( observeList[t](j,c) ) ) {
            obs_t += 1;
            c_cov += matrixList[t](i,c) * matrixList[t](j,c);
          }
          
        }
        
        if (obs_t >= 1) {
          cov_acc += (c_cov/obs_t);
        }
      }
      
      weight_cov(i,j) = cov_acc;
    }
  }
  
  
  // Use SelfAdjointView to copy the upper triangular part to the lower triangular part
  //weight_cov.triangularView<Eigen::Lower>() = weight_cov.transpose().selfadjointView<Eigen::Upper>();
  Eigen::MatrixXd weight_cov_t = weight_cov.transpose();
  weight_cov += weight_cov_t;
  
  for (int i = 0; i < dimkm1; i++){
    weight_cov(i,i) /= 2;
  }
  
  //for debug purpose
  //double ele = k_mat(0,0);
  //bool ele = std::isnan( tensor(0,0,0,0) );
  //return ele;
  //return matrixList[1];
  
  return weight_cov;
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//



/*** R
#K1_cov_est(1:24, c(2,12), 2)
#should be
#[,1] [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#[1,]  2.5  5.5   8.5  11.5  14.5  17.5  20.5  23.5  26.5  29.5  32.5  35.5
#[2,]  5.5 12.5  19.5  26.5  33.5  40.5  47.5  54.5  61.5  68.5  75.5  82.5
#[3,]  8.5 19.5  30.5  41.5  52.5  63.5  74.5  85.5  96.5 107.5 118.5 129.5
#[4,] 11.5 26.5  41.5  56.5  71.5  86.5 101.5 116.5 131.5 146.5 161.5 176.5
#[5,] 14.5 33.5  52.5  71.5  90.5 109.5 128.5 147.5 166.5 185.5 204.5 223.5
#[6,] 17.5 40.5  63.5  86.5 109.5 132.5 155.5 178.5 201.5 224.5 247.5 270.5
#[7,] 20.5 47.5  74.5 101.5 128.5 155.5 182.5 209.5 236.5 263.5 290.5 317.5
#[8,] 23.5 54.5  85.5 116.5 147.5 178.5 209.5 240.5 271.5 302.5 333.5 364.5
#[9,] 26.5 61.5  96.5 131.5 166.5 201.5 236.5 271.5 306.5 341.5 376.5 411.5
#[10,] 29.5 68.5 107.5 146.5 185.5 224.5 263.5 302.5 341.5 380.5 419.5 458.5
#[11,] 32.5 75.5 118.5 161.5 204.5 247.5 290.5 333.5 376.5 419.5 462.5 505.5
#[12,] 35.5 82.5 129.5 176.5 223.5 270.5 317.5 364.5 411.5 458.5 505.5 552.5

*/
