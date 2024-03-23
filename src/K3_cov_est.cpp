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
Eigen::MatrixXd K3_cov_est(NumericVector x, NumericVector dim, const int k) {
  // NumericVector x = clone(xi); (for debug, no need now)
  // k is mode starting from 2 (using the way that R counts)
  // since time mode is mode 1;
  // so for K=4, we have k=2,3,4 ;

  int dim1 = dim[0];
  int dim2 = dim[1];
  int dim3 = dim[2];
  int dim4 = dim[3];
  int dimkm1 = dim[k-1];
  
  Eigen::Tensor<double, 4> tensor(dim1, dim2, dim3, dim4);
  
  for (int i=0; i < dim1; i++){
    for (int j=0; j < dim2; j++){
      for (int s=0; s < dim3; s++){
        for (int l=0; l < dim4; l++){
          tensor(i,j,s,l) = x(i + j*dim1 + s*dim1*dim2 + l*dim1*dim2*dim3);
        }
      }
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
        //Eigen::DSizes<int, 4> indices;
        Eigen::array<Eigen::Index, 4> indices;
        indices[0] = t;
        indices[k-1] = u1;
        for (int w = 1; w < 4; w++) {
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
#K3_cov_est(c(1:478, NA, NA), c(20,3,4,2), 2)
#should be
# [,1]     [,2]     [,3]
#[1,] 540428.0 575708.0 610081.7
#[2,] 575708.0 614188.0 651741.7
#[3,] 610081.7 651741.7 693401.7
#R runs in: 0.495352 secs

#K3_cov_est(c(1:478, NA, NA), c(20,3,4,2), 3)
# should be
# [,1]     [,2]     [,3]     [,4]
#[1,] 224101.0 278281.0 332461.0 385874.7
#[2,] 278281.0 354061.0 429841.0 504794.7
#[3,] 332461.0 429841.0 527221.0 623714.7
#[4,] 385874.7 504794.7 623714.7 742634.7

#K3_cov_est(c(1:478, NA, NA), c(20,3,4,2), 4)
# should be
#[,1]      [,2]
#[1,] 231842.0  578175.7
#[2,] 578175.7 1616175.7


#dt <- array(c(12,NA,10,NA,1:5,NA, 23:32, NA) , c(2,3,2,2))
#K3_cov_est(c(dt), dim(dt), 2)
# [,1]   [,2]   [,3]
#[1,] 1799.5 1220.5 1167.0
#[2,] 1220.5 1025.5 1056.5
#[3,] 1167.0 1056.5 1525.5

#K3_cov_est(c(dt), dim(dt), 3)
#[,1]   [,2]
#[1,] 2524 1561.0
#[2,] 1561 1826.5

#K3_cov_est(c(dt), dim(dt), 4)
# [,1] [,2]
#[1,] 836.5  965
#[2,] 965.0 3514
*/
