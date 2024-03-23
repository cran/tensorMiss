#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <math.h>
#include <vector>
#include <assert.h>

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

void sortrows(std::vector<std::vector<double>>& matrix, int col) {    
  std::sort(matrix.begin(),
            matrix.end(),
            [col](const std::vector<double>& lhs, const std::vector<double>& rhs) {
              return lhs[col] > rhs[col];
            });
}

// template for quantile (assuming sorted in descending order)
template <typename t1, typename t2> typename t1::value_type quant(const t1 &x, t2 q)
{
  assert(q >= 0.0 && q <= 1.0);
  
  const auto n  = x.size();
  const auto id = (n - 1) * (1-q);
  const auto lo = floor(id);
  const auto hi = ceil(id);
  const auto qs = x[lo];
  const auto h  = (id - lo);
  
  return (1.0 - h) * qs + h * x[hi];
}





// [[Rcpp::export]]
double partition_MSE(std::vector<double> x1,
                     std::vector<double> x2, int par_num = 100) {
  // x1: observed vector of length of n
  // x2: imputed vector of length of n
  
  int N = x1.size();
  
  // declare sortrows function
  void sortrows(std::vector<std::vector<double>>& matrix, int col);
  
  std::vector<std::vector<double>> dt_mat(N);
  for (int i = 0; i < N; i++){
    std::vector<double> dt_i(2);
    dt_i[0] = x1[i];
    dt_i[1] = x2[i];
    dt_mat[i] = dt_i;
  }

  sortrows(dt_mat, 0);

  
  // extract observed vector and imputed vector after sorted
  std::vector<float> quan_x_i1(N);
  for (int i = 0; i < N; i++){
    quan_x_i1[i] = dt_mat[i][0];
  }
  std::vector<float> quan_x_i2(N);
  for (int i = 0; i < N; i++){
    quan_x_i2[i] = dt_mat[i][1];
  }
  
  // obtain quantile vector
  std::vector<float> quantile_vec(par_num);
  for (int i = 0; i < par_num; i++){
    float quan_i;
    quan_i = (float)i / (float)par_num;
    quantile_vec[i] = quant(quan_x_i1, quan_i);
  }
  
  // for debug
  // for (int i=0; i< quantile_vec.size(); i++){
  //   std::cout << quantile_vec[i] << std::endl;
  // }
  
  double y1 = 0;
  double y2 = 0;
  
  int rec1 = 0;
  for (int i = par_num - 1; i >= 0; i--){
    double y_i = 0;
    
    for (int j = rec1; j < N; j++){
      if ( quan_x_i1[j] >= quantile_vec[i] ){
        y_i += quan_x_i1[j] - quan_x_i2[j];
      } else {
        rec1 = j;
        break;
      }
    }
    
    y1 += pow(y_i, 2);
  }
  
  int rec2 = 0;
  for (int i = par_num - 1; i >= 0; i--){
    double y_i = 0;
    
    for (int j = rec2; j < N; j++){
      if ( quan_x_i1[j] >= quantile_vec[i] ){
        y_i += quan_x_i1[j];
      } else {
        rec2 = j;
        break;
      }
    }
    
    y2 += pow(y_i, 2);
  }
  
  
  double result;
  result = y1 / y2;

  
  return result;
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//



/*** R
# dt_1 <- c(2, 3, 7, 1)
# dt_2 <- c(-2, 0.5, -8, 9)
# partition_MSE(dt_1, dt_2, 1)
# should be: (sum(dt_1-dt_2)^2)/(sum(dt_1)^2) = 1.078

# dt_1 <- c(2, 3, 4, 1)
# dt_2 <- c(2, 3, 4, 1)
# partition_MSE(dt_1, dt_2, 4)
# should be: 0

# start_time <- Sys.time()
# partition_MSE(c(ten_comp), c(Y.imp))
# end_time <- Sys.time()
# print(end_time-start_time)

*/
