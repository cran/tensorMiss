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
Eigen::VectorXd K3_Ft_est(NumericVector x, NumericMatrix A1, NumericMatrix A2, NumericMatrix A3) {
  // NumericVector x = clone(xi); (for debug, no need now)
  
  // declare my_kroneckerProduct function
  Eigen::MatrixXd my_kroneckerProduct(NumericMatrix A, NumericMatrix B);
  
  // design matrix
  Eigen::MatrixXd kron = my_kroneckerProduct(A3, A2);
  kron = my_kroneckerProduct(wrap(kron), A1);
  
  // record NA
  LogicalVector miss_bool = !is_na(x);
  LogicalVector col_bool(kron.cols(), true);
  
  
  Eigen::VectorXd Y_reg = as<Eigen::VectorXd>(x[miss_bool]);
  
  // Create a new Eigen matrix containing only the selected rows
  long num_true = std::count(miss_bool.begin(), miss_bool.end(), true);
  Eigen::MatrixXd kron_cut(num_true, kron.cols());
  NumericVector true_index(num_true);
  
  long j = 0;
  for (long i = 0; i < miss_bool.size(); i++) {
    if (miss_bool(i)) {
      //true_index.push_back(i);
      true_index(j) = i;
      j += 1;
    }
  }
  
  for (long i = 0; i < num_true; i++) {
    kron_cut.row(i) = kron.row( true_index(i) );
  }
  
  
  Eigen::MatrixXd kron_inv = (kron_cut.transpose() * kron_cut).inverse() * kron_cut.transpose();
  Eigen::VectorXd result_reg = kron_inv * Y_reg;
  
  
  Eigen::VectorXd result = result_reg;
  return result;
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//



/*** R

# isTRUE(all.equal( F.est[1:length(c(F.est))], F.est_compare[1:length(c(F.est_compare))] ))

*/
