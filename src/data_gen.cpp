#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <math.h>
#include <cmath>
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

Eigen::MatrixXd my_kroneckerProduct(NumericMatrix A, NumericMatrix B) {
  // Convert NumericMatrix to Eigen matrices
  Eigen::Map<Eigen::MatrixXd> A_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(A));
  Eigen::Map<Eigen::MatrixXd> B_eigen(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(B));
  
  // Calculate the Kronecker product using Eigen
  Eigen::MatrixXd result;
  result = Eigen::kroneckerProduct(A_eigen, B_eigen);
  
  return result;
}





Eigen::Tensor<double, 2> core_1gen(int n, IntegerVector r,
                                   NumericVector coef,
                                   std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // K: order of tensor at time t (we only write up to K=4), here is K=1, so omit parameter K;
  // n: total number of time points;
  // r: rank of core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(r(0));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  Eigen::Tensor<double, 2> tensor(n, r(0));
  
  for (int i=0; i < r(0); i++){
    // ARIMA AR coefficients
    Eigen::VectorXd ar_coeffs(5);
    for (int i=0; i < 5; i++){
      if ( (i+1) > coef.size()){
        ar_coeffs[i] = 0;
      } else {
        ar_coeffs[i] = coef[i]; 
      }
    }
    
    Eigen::VectorXd ar_dt(n);
    ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
    
    rec += 1;
    
    for (int t=0; t < n; t++){
      tensor(t,i) = ar_dt(t);
    }
  }

  //auto tensorx = tensor.reshape(tensor.size());
  return tensor;
}

Eigen::Tensor<double, 3> core_2gen(int n, IntegerVector r,
                                   NumericVector coef,
                                   std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // K: order of tensor at time t (we only write up to K=4), here is K=2, so omit parameter K;
  // n: total number of time points;
  // r: rank of core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(r(0)*r(1));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  Eigen::Tensor<double, 3> tensor(n, r(0), r(1));
  
  for (int i=0; i < r(0); i++){
    for (int j=0; j < r(1); j++){
      // ARIMA AR coefficients
      Eigen::VectorXd ar_coeffs(5);
      for (int i=0; i < 5; i++){
        if ( (i+1) > coef.size()){
          ar_coeffs[i] = 0;
        } else {
          ar_coeffs[i] = coef[i]; 
        }
      }
      
      Eigen::VectorXd ar_dt(n);
      ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
      
      rec += 1;
      
      for (int t=0; t < n; t++){
        tensor(t,i,j) = ar_dt(t);
      }
    }
  }
  
  //auto tensorx = tensor.reshape(tensor.size());
  return tensor;
}


Eigen::Tensor<double, 4> core_3gen(int n, IntegerVector r,
                                   NumericVector coef,
                                   std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // K: order of tensor at time t (we only write up to K=4), here is K=3, so omit parameter K;
  // n: total number of time points;
  // r: rank of core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(r(0)*r(1)*r(2));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  Eigen::Tensor<double, 4> tensor(n, r(0), r(1), r(2));
  
  for (int i=0; i < r(0); i++){
    for (int j=0; j < r(1); j++){
      for (int a=0; a < r(2); a++){
        // ARIMA AR coefficients
        Eigen::VectorXd ar_coeffs(5);
        for (int i=0; i < 5; i++){
          if ( (i+1) > coef.size()){
            ar_coeffs[i] = 0;
          } else {
            ar_coeffs[i] = coef[i]; 
          }
        }
        
        Eigen::VectorXd ar_dt(n);
        ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);

        rec += 1;
        
        for (int t=0; t < n; t++){
          tensor(t,i,j,a) = ar_dt(t);
        }
      }
    }
  }
  
  //auto tensorx = tensor.reshape(tensor.size());
  return tensor;
}



Eigen::Tensor<double, 5> core_4gen(int n, IntegerVector r,
                                   NumericVector coef,
                                   std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // K: order of tensor at time t (we only write up to K=4), here is K=4, so omit parameter K;
  // n: total number of time points;
  // r: rank of core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(r(0)*r(1)*r(2)*r(3));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  Eigen::Tensor<double, 5> tensor(n, r(0), r(1), r(2), r(3));
  
  for (int i=0; i < r(0); i++){
    for (int j=0; j < r(1); j++){
      for (int a=0; a < r(2); a++){
        for (int b=0; b < r(3); b++){
          // ARIMA AR coefficients
          Eigen::VectorXd ar_coeffs(5);
          for (int i=0; i < 5; i++){
            if ( (i+1) > coef.size()){
              ar_coeffs[i] = 0;
            } else {
              ar_coeffs[i] = coef[i]; 
            }
          }
          
          Eigen::VectorXd ar_dt(n);
          ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
          
          rec += 1;
          
          for (int t=0; t < n; t++){
            tensor(t,i,j,a,b) = ar_dt(t);
          } 
        }
      }
    }
  }
  //auto tensorx = tensor.reshape(tensor.size());
  return tensor;
}



Eigen::Tensor<double, 2> common_e_1gen(int n, IntegerVector re,
                                       NumericVector coef,
                                       std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // re: rank of error core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(re(0));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  Eigen::Tensor<double, 2> tensor(n, re(0));
  
  for (int i=0; i < re(0); i++){
    // ARIMA AR coefficients
    Eigen::VectorXd ar_coeffs(5);
    for (int i=0; i < 5; i++){
      if ( (i+1) > coef.size()){
        ar_coeffs[i] = 0;
      } else {
        ar_coeffs[i] = coef[i]; 
      }
    }
    
    Eigen::VectorXd ar_dt(n);
    ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
    
    rec += 1;
    
    for (int t=0; t < n; t++){
      tensor(t,i) = ar_dt(t);
    }
  }
  
  //auto tensorx = tensor.reshape(tensor.size());
  return tensor;
}

Eigen::Tensor<double, 3> common_e_2gen(int n, IntegerVector re,
                                       NumericVector coef,
                                       std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // re: rank of error core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(re(0)*re(1));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  Eigen::Tensor<double, 3> tensor(n, re(0), re(1));
  
  for (int i=0; i < re(0); i++){
    for (int j=0; j < re(1); j++){
      // ARIMA AR coefficients
      Eigen::VectorXd ar_coeffs(5);
      for (int i=0; i < 5; i++){
        if ( (i+1) > coef.size()){
          ar_coeffs[i] = 0;
        } else {
          ar_coeffs[i] = coef[i]; 
        }
      }
      
      Eigen::VectorXd ar_dt(n);
      ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
      
      rec += 1;
      
      for (int t=0; t < n; t++){
        tensor(t,i,j) = ar_dt(t);
      }
    }
  }
  
  //auto tensorx = tensor.reshape(tensor.size());
  return tensor;
}


Eigen::Tensor<double, 4> common_e_3gen(int n, IntegerVector re,
                                       NumericVector coef,
                                       std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // re: rank of error core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(re(0)*re(1)*re(2));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  Eigen::Tensor<double, 4> tensor(n, re(0), re(1), re(2));
  
  for (int i=0; i < re(0); i++){
    for (int j=0; j < re(1); j++){
      for (int a=0; a < re(2); a++){
        // ARIMA AR coefficients
        Eigen::VectorXd ar_coeffs(5);
        for (int i=0; i < 5; i++){
          if ( (i+1) > coef.size()){
            ar_coeffs[i] = 0;
          } else {
            ar_coeffs[i] = coef[i]; 
          }
        }
        
        Eigen::VectorXd ar_dt(n);
        ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
        
        rec += 1;
        
        for (int t=0; t < n; t++){
          tensor(t,i,j,a) = ar_dt(t);
        }
      }
    }
  }
  
  //auto tensorx = tensor.reshape(tensor.size());
  return tensor;
}



Eigen::Tensor<double, 5> common_e_4gen(int n, IntegerVector re,
                                       NumericVector coef,
                                       std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // re: rank of error core tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(re(0)*re(1)*re(2)*re(3));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  
  Eigen::Tensor<double, 5> tensor(n, re(0), re(1), re(2), re(3));
  
  for (int i=0; i < re(0); i++){
    for (int j=0; j < re(1); j++){
      for (int a=0; a < re(2); a++){
        for (int b=0; b < re(3); b++){
          // ARIMA AR coefficients
          Eigen::VectorXd ar_coeffs(5);
          for (int i=0; i < 5; i++){
            if ( (i+1) > coef.size()){
              ar_coeffs[i] = 0;
            } else {
              ar_coeffs[i] = coef[i]; 
            }
          }
          
          Eigen::VectorXd ar_dt(n);
          ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
          
          rec += 1;
          
          for (int t=0; t < n; t++){
            tensor(t,i,j,a,b) = ar_dt(t);
          } 
        }
      }
    }
  }
  return tensor;
}




Eigen::Tensor<double, 2> idio_e_1gen(int n, IntegerVector d,
                                     NumericVector coef,
                                     std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // d: dimension of data tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(d(0));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  
  //std::random_device rd;
  //std::mt19937 generator(rd());
  std::normal_distribution<double> gauss_dist(0.0, 1.0);
  
  Eigen::Tensor<double, 2> tensor(n, d(0));
  
  for (int i=0; i < d(0); i++){
    // ARIMA AR coefficients
    Eigen::VectorXd ar_coeffs(5);
    for (int i=0; i < 5; i++){
      if ( (i+1) > coef.size()){
        ar_coeffs[i] = 0;
      } else {
        ar_coeffs[i] = coef[i]; 
      }
    }
    
    Eigen::VectorXd ar_dt(n);
    ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
    
    rec += 1;
    
    for (int t=0; t < n; t++){
      // tensor(t,i) = ar_dt(t) * abs(gauss_dist(generator));
      tensor(t,i) = ar_dt(t);
    }
  }
  
  return tensor;
}

Eigen::Tensor<double, 3> idio_e_2gen(int n, IntegerVector d,
                                     NumericVector coef,
                                     std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // d: dimension of data tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(d(0)*d(1));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  
  //std::random_device rd;
  //std::mt19937 generator(rd());
  std::normal_distribution<double> gauss_dist(0.0, 1.0);
  
  Eigen::Tensor<double, 3> tensor(n, d(0), d(1));
  
  for (int i=0; i < d(0); i++){
    for (int j=0; j < d(1); j++){
      // ARIMA AR coefficients
      Eigen::VectorXd ar_coeffs(5);
      for (int i=0; i < 5; i++){
        if ( (i+1) > coef.size()){
          ar_coeffs[i] = 0;
        } else {
          ar_coeffs[i] = coef[i]; 
        }
      }
      
      Eigen::VectorXd ar_dt(n);
      ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
      
      rec += 1;
      
      for (int t=0; t < n; t++){
        // tensor(t,i,j) = ar_dt(t) * abs(gauss_dist(generator));
        tensor(t,i,j) = ar_dt(t);
      }
    }
  }
  
  return tensor;
}


Eigen::Tensor<double, 4> idio_e_3gen(int n, IntegerVector d,
                                     NumericVector coef,
                                     std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // d: dimension of data tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(d(0)*d(1)*d(2));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  
  //std::random_device rd;
  //std::mt19937 generator(rd());
  std::normal_distribution<double> gauss_dist(0.0, 1.0);
  
  Eigen::Tensor<double, 4> tensor(n, d(0), d(1), d(2));
  
  for (int i=0; i < d(0); i++){
    for (int j=0; j < d(1); j++){
      for (int a=0; a < d(2); a++){
        // ARIMA AR coefficients
        Eigen::VectorXd ar_coeffs(5);
        for (int i=0; i < 5; i++){
          if ( (i+1) > coef.size()){
            ar_coeffs[i] = 0;
          } else {
            ar_coeffs[i] = coef[i]; 
          }
        }
        
        Eigen::VectorXd ar_dt(n);
        ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
        
        rec += 1;
        
        for (int t=0; t < n; t++){
          // tensor(t,i,j,a) = ar_dt(t) * abs(gauss_dist(generator));
          tensor(t,i,j,a) = ar_dt(t);
        }
      }
    }
  }
  
  return tensor;
}



Eigen::Tensor<double, 5> idio_e_4gen(int n, IntegerVector d,
                                     NumericVector coef,
                                     std::mt19937 generator, bool heavy_tailed = false, int t_df = 5){
  // n: total number of time points;
  // d: dimension of data tensors, written in a vector of length K;
  
  // declare ar_sim function
  Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                         bool heavy_tailed = false, int t_df = 5);
  
  std::vector<std::mt19937> genList(d(0)*d(1)*d(2)*d(3));
  std::uniform_int_distribution<int>  distr(1, 100000000);
  for (int i=0; i < genList.size(); i++){
    std::mt19937 gen_i(distr(generator));
    genList[i] = gen_i;
  }
  int rec = 0;
  
  //std::random_device rd;
  //std::mt19937 generator(rd());
  std::normal_distribution<double> gauss_dist(0.0, 1.0);
  
  Eigen::Tensor<double, 5> tensor(n, d(0), d(1), d(2), d(3));
  
  for (int i=0; i < d(0); i++){
    for (int j=0; j < d(1); j++){
      for (int a=0; a < d(2); a++){
        for (int b=0; b < d(3); b++){
          // ARIMA AR coefficients
          Eigen::VectorXd ar_coeffs(5);
          for (int i=0; i < 5; i++){
            if ( (i+1) > coef.size()){
              ar_coeffs[i] = 0;
            } else {
              ar_coeffs[i] = coef[i]; 
            }
          }
          
          Eigen::VectorXd ar_dt(n);
          ar_dt = ar_sim(n, ar_coeffs, genList[rec], heavy_tailed, t_df);
          
          rec += 1;
          
          for (int t=0; t < n; t++){
            // tensor(t,i,j,a,b) = ar_dt(t) * abs(gauss_dist(generator));
            tensor(t,i,j,a,b) = ar_dt(t);
          }
        }
      }
    }
  }
  return tensor;
}




Eigen::VectorXd ar_sim(int n, Eigen::VectorXd ar_coeffs, std::mt19937 generator,
                       bool heavy_tailed = false, int t_df = 5){
  
  int n_start = ar_coeffs.size();
  int p = ar_coeffs.size();
  
  // Simulate innovations
  //std::random_device rd;  // Obtain a random seed from the hardware
  //std::mt19937 generator(rd()); // Mersenne Twister pseudo-random generator
  
  Eigen::VectorXd innovations(n+n_start);
  
  if (heavy_tailed) {
    std::student_t_distribution<double> t_dist(t_df);
    for (int i = 0; i < innovations.size(); i++) {
      innovations[i] = t_dist(generator);
    }
  } else {
    std::normal_distribution<double> gauss_dist(0.0, 1.0);
    // std::normal_distribution<double> gauss_dist(0.0, 0.5);
    for (int i = 0; i < innovations.size(); i++) {
      innovations[i] = gauss_dist(generator);
    }
  }
  
  // Initialize the time series
  Eigen::VectorXd seriesx(n+n_start);
  for (int i = 0; i < n_start; i++) {
    seriesx(i) = 0.0;
  }
  
  // Generate the time series data
  for (int i = n_start; i < seriesx.size(); i++) {
    double value = 0.0;
    
    for (int j = 0; j < p; j++) {
      value += ar_coeffs(j) * seriesx(i - j - 1);
    }
    
    value += innovations(i);
    seriesx(i) = value;
  }
  
  Eigen::VectorXd series(n);
  for (int i = 0; i < n; i++){
    series[i] = seriesx(i+n_start);
  }
  
  // standardise the series
  double sample_var = 0;
  for (int i = 0; i < n; i++){
    sample_var += pow(series[i], 2);
  }
  
  sample_var /= n;
  series /= pow(sample_var, 0.5);
  
  return series;
}








// [[Rcpp::export]]
List data_gen(int K, //   the number of modes for the tensor time series
              int n, //   length of time series
              IntegerVector d, //   dimensions in each mode, written in a vector of length K. e.g. c(40,40)
              IntegerVector r, //   rank of core tensors, written in a vector of length K. e.g. c(2,2)
              IntegerVector re, //  rank of the cross-sectional common error core tensors, written in a vector of length K. e.g. c(2,2)
              List eta, // quantities controlling factor strengths in A_k, written in a list of K vectors. e.g. eta = list(c(0,0.2),c(0,0.2))
              NumericVector coef_f, //   quantities controlling AR(5) coefficients in F_t (\phi_1, \phi_2, \phi_3, \phi_4, \phi_5), written in a vector of length 5. e.g. coef_f = c(0.7, 0.3, -0.4, 0.2, -0.1)
              NumericVector coef_fe, //   quantities controlling AR(5) coefficients in F_{e,t} (\phi_1, \phi_2, \phi_3, \phi_4, \phi_5), written in a vector of length 5. e.g. coef_fe = c(-0.7, -0.3, -0.4, 0.2, 0.1)
              NumericVector coef_e, //   quantities controlling AR(5) coefficients in \epsilon_t (\phi_1, \phi_2, \phi_3, \phi_4, \phi_5), written in a vector of length 5. e.g. coef_e = c(0.8, 0.4, -0.4, 0.2, -0.1)
              bool heavy_tailed = false, // whether to generate data from heavy tailed distribution. If FALSE, generate from N(0,1); if TRUE, generate from t-distribution
              int t_df = 3, // if heavy_tailed = TRUE, the degree of freedom for t-distribution
              int seed = 2023 //set random seed
) {
  // output : a list of six elements:
  // output$A: a list of K factor loading matrices
  // output$Ao: a matrix of kronecker product of the K loading matrices, written as A_K \otimes ... \otimes A_1
  // output$F_ts: vectorisation of time series of the core tensors
  // output$E_common: vectorisation of the error common component core tensors
  // output$Aeo: a matrix of kronecker product of the K error loading matrices, written as A_{e,K} \otimes ... \otimes A_{e,1}
  // output$E_idio: vectorisation of the error idiosyncratic tensors
  
  
  std::mt19937 generator1(seed);
  
  
  // Generate the core tensor (independent AR process) ------------------------------
  // declare core_gen functions
  Eigen::Tensor<double, 2> core_1gen(int n, IntegerVector r, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 3> core_2gen(int n, IntegerVector r, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 4> core_3gen(int n, IntegerVector r, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 5> core_4gen(int n, IntegerVector r, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  
  Eigen::VectorXd f_series;
  
  if (K == 1){
    Eigen::Tensor<double, 2> f_seriesx(n, r(0));
    f_seriesx = core_1gen(n, r, coef_f, generator1, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> f_series_pre(f_seriesx.data(), f_seriesx.size());
    f_series = f_series_pre;
  } else if (K == 2){
    Eigen::Tensor<double, 3> f_seriesx(n, r(0), r(1));
    f_seriesx = core_2gen(n, r, coef_f, generator1, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> f_series_pre(f_seriesx.data(), f_seriesx.size());
    f_series = f_series_pre;
  } else if (K == 3){
    Eigen::Tensor<double, 4> f_seriesx(n, r(0), r(1), r(2));
    f_seriesx = core_3gen(n, r, coef_f, generator1, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> f_series_pre(f_seriesx.data(), f_seriesx.size());
    f_series = f_series_pre;
  } else if (K == 4){
    Eigen::Tensor<double, 5> f_seriesx(n, r(0), r(1), r(2), r(3));
    f_seriesx = core_4gen(n, r, coef_f, generator1, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> f_series_pre(f_seriesx.data(), f_seriesx.size());
    f_series = f_series_pre;
  }
  

  
  // Generate factor loading matrices ------------------------------------------------------------
  std::mt19937 generator(seed-1);
  
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> A(K);
  
  //std::default_random_engine generator;
  // std::random_device rd;  // Obtain a random seed from the hardware
  // std::mt19937 generator(rd()); // Mersenne Twister pseudo-random generator
    
    for (int k = 0; k < K; k++)
    {
      std::normal_distribution<double> gauss_dist(0.0, 1.0);

      Eigen::VectorXd U_k_vec( d(k)*r(k) );
      for (int i = 0; i < U_k_vec.size(); i++) {
        //U_k_vec[i] = unif_dist(generator);
        U_k_vec[i] = gauss_dist(generator);
      }

      Eigen::Map<Eigen::MatrixXd> U_k(U_k_vec.data(), d(k), r(k));

      Eigen::VectorXd eta_k = eta(k);
      
      Eigen::VectorXd R_k_vec(r(k));
      for (int i = 0; i < r(k); i++){
        R_k_vec[i] = pow(d(k), -eta_k(i));
      }
      
      //Eigen::DiagonalMatrix<double, r(k)> R_k(R_k_vec);
      Eigen::MatrixXd R_k = Eigen::MatrixXd::Constant(r(k), r(k), 0.0);
      for (int i = 0; i < r(k); i++){
        R_k(i,i) = R_k_vec(i);
      }
      
      Eigen::MatrixXd A_k = U_k * R_k;
      A[k] = A_k;
    }
  
  // Generate common components ------------------------------------------------------------
  // declare my_kroneckerProduct function
  Eigen::MatrixXd my_kroneckerProduct(NumericMatrix A, NumericMatrix B);
  
  Eigen::MatrixXd A_otimes;
  if (K == 1){
    A_otimes = A[0];
  } else if (K == 2){
    A_otimes = A[0];
    A_otimes = my_kroneckerProduct(wrap(A[1]), wrap(A_otimes));
  } else if (K == 3){
    A_otimes = A[0];
    A_otimes = my_kroneckerProduct(wrap(A[1]), wrap(A_otimes));
    A_otimes = my_kroneckerProduct(wrap(A[2]), wrap(A_otimes));
  } else if (K == 4){
    A_otimes = A[0];
    A_otimes = my_kroneckerProduct(wrap(A[1]), wrap(A_otimes));
    A_otimes = my_kroneckerProduct(wrap(A[2]), wrap(A_otimes));
    A_otimes = my_kroneckerProduct(wrap(A[3]), wrap(A_otimes));
  }
  
  
  // Generate the error tensors ------------------------------------------------------------
  std::mt19937 generator2(seed+1);
  // declare common_e_gen functions
  Eigen::Tensor<double, 2> common_e_1gen(int n, IntegerVector re, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 3> common_e_2gen(int n, IntegerVector re, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 4> common_e_3gen(int n, IntegerVector re, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 5> common_e_4gen(int n, IntegerVector re, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  
  Eigen::VectorXd common_e;
  
  if (K == 1){
    Eigen::Tensor<double, 2> common_e_x(n, re(0));
    common_e_x = common_e_1gen(n, re, coef_fe, generator2, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> common_e_pre(common_e_x.data(), common_e_x.size());
    common_e = common_e_pre;
  } else if (K == 2){
    Eigen::Tensor<double, 3> common_e_x(n, re(0), re(1));
    common_e_x = common_e_2gen(n, re, coef_fe, generator2, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> common_e_pre(common_e_x.data(), common_e_x.size());
    common_e = common_e_pre;
  } else if (K == 3){
    Eigen::Tensor<double, 4> common_e_x(n, re(0), re(1), re(2));
    common_e_x = common_e_3gen(n, re, coef_fe, generator2, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> common_e_pre(common_e_x.data(), common_e_x.size());
    common_e = common_e_pre;
  } else if (K == 4){
    Eigen::Tensor<double, 5> common_e_x(n, re(0), re(1), re(2), re(3));
    common_e_x = common_e_4gen(n, re, coef_fe, generator2, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> common_e_pre(common_e_x.data(), common_e_x.size());
    common_e = common_e_pre;
  }
  
  
  // Generate the error loading matrices ------------------------------------------------------------
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> Ae(K);
  
  //std::default_random_engine generator;
  //std::random_device rd_e;  // Obtain a random seed from the hardware
  //std::mt19937 generator_e(rd_e()); // Mersenne Twister pseudo-random generator
  //std::random_device rd_sparse_e;  // Obtain a random seed from the hardware
  //std::mt19937 generator_sparse_e(rd_sparse_e()); // Mersenne Twister pseudo-random generator
  
  for (int k = 0; k < K; k++)
  {
    std::normal_distribution<double> gauss_dist(0.0, 1.0);
    std::bernoulli_distribution bernoulli_dist(0.95);
    
    std::vector<bool> Ae_k_sparse_vec( d(k) * re(k) );
    for (int i = 0; i < Ae_k_sparse_vec.size(); i++) {
      Ae_k_sparse_vec[i] = bernoulli_dist(generator1);
      //Ae_k_sparse_vec[i] = bernoulli_dist(generator_sparse_e);
    }
    
    Eigen::VectorXd Ae_k_vec( d(k) * re(k) );
    for (int i = 0; i < Ae_k_vec.size(); i++) {
      if ( Ae_k_sparse_vec[i] ){
        Ae_k_vec[i] = 0;
      } else {
        Ae_k_vec[i] = gauss_dist(generator1);
        //Ae_k_vec[i] = gauss_dist(generator_e);
      }
    }
    
    Eigen::Map<Eigen::MatrixXd> Ae_k(Ae_k_vec.data(), d(k), re(k));
    
    Ae[k] = Ae_k;
  }
  
  Eigen::MatrixXd Ae_otimes;
  if (K == 1){
    Ae_otimes = Ae[0];
  } else if (K == 2){
    Ae_otimes = Ae[0];
    Ae_otimes = my_kroneckerProduct(wrap(Ae[1]), wrap(Ae_otimes));
  } else if (K == 3){
    Ae_otimes = Ae[0];
    Ae_otimes = my_kroneckerProduct(wrap(Ae[1]), wrap(Ae_otimes));
    Ae_otimes = my_kroneckerProduct(wrap(Ae[2]), wrap(Ae_otimes));
  } else if (K == 4){
    Ae_otimes = Ae[0];
    Ae_otimes = my_kroneckerProduct(wrap(Ae[1]), wrap(Ae_otimes));
    Ae_otimes = my_kroneckerProduct(wrap(Ae[2]), wrap(Ae_otimes));
    Ae_otimes = my_kroneckerProduct(wrap(Ae[3]), wrap(Ae_otimes));
  }
  
  
  // Generate idiosyncratic error tensors ------------------------------------------------------------
  std::mt19937 generator3(seed+2);
  // declare common_e_gen functions
  Eigen::Tensor<double, 2> idio_e_1gen(int n, IntegerVector d, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 3> idio_e_2gen(int n, IntegerVector d, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 4> idio_e_3gen(int n, IntegerVector d, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  Eigen::Tensor<double, 5> idio_e_4gen(int n, IntegerVector d, NumericVector coef, std::mt19937 generator, bool heavy_tailed = false, int t_df = 5);
  
  Eigen::VectorXd idio_e;
  
  if (K == 1){
    Eigen::Tensor<double, 2> idio_e_x(n, d(0));
    idio_e_x = idio_e_1gen(n, d, coef_e, generator3, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> idio_e_pre(idio_e_x.data(), idio_e_x.size());
    idio_e = idio_e_pre;
  } else if (K == 2){
    Eigen::Tensor<double, 3> idio_e_x(n, d(0), d(1));
    idio_e_x = idio_e_2gen(n, d, coef_e, generator3, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> idio_e_pre(idio_e_x.data(), idio_e_x.size());
    idio_e = idio_e_pre;
  } else if (K == 3){
    Eigen::Tensor<double, 4> idio_e_x(n, d(0), d(1), d(2));
    idio_e_x = idio_e_3gen(n, d, coef_e, generator3, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> idio_e_pre(idio_e_x.data(), idio_e_x.size());
    idio_e = idio_e_pre;
  } else if (K == 4){
    Eigen::Tensor<double, 5> idio_e_x(n, d(0), d(1), d(2), d(3));
    idio_e_x = idio_e_4gen(n, d, coef_e, generator3, heavy_tailed, t_df);
    Eigen::Map<Eigen::VectorXd> idio_e_pre(idio_e_x.data(), idio_e_x.size());
    idio_e = idio_e_pre;
  }
  
  
  
  
  
  
  
  // return List::create(Named("A")= A);
  return List::create(Named("A")= A,
                      Named("Ao")= A_otimes,
                      Named("F_ts")= f_series,
                      Named("E_common")= common_e,
                      Named("Aeo")= Ae_otimes,
                      Named("E_idio")= idio_e);
  }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//



/*** R
# K = 3
# 
# I1 = 30
# I2 = 30
# I3 = 30
# r1 = 2
# r2 = 3
# r3 = 3
# 
# TT = 20
# d = c(I1,I2,I3)
# r = c(r1,r2,r3)
# re = c(2,2,2)
# eta = list(c(0.3,0.0), c(0.0,0.0, 0.0), c(0.0,0.0, 0.0))
# coef_f = c(0.7, 0.3, -0.4, 0.2, -0.1)
# coef_fe = c(-0.7, -0.3, -0.4, 0.2, 0.1)
# coef_e = c(0.8, 0.4, -0.4, 0.2, -0.1)
# # coef_f = numeric()
# # coef_fe = numeric() 
# # coef_e = numeric()
# #set.seed(10)
# data_test = data_gen(K,TT,d,r,re,eta, coef_f, coef_fe, coef_e, heavy_tailed=TRUE, t_df=6, seed = 2023)
*/
