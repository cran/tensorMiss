
#' HAC covariance estimator for asymptotic normality on each row j of loading matrix estimator
#'
#' @description Computing the HAC covariance estimator for asymptotic normality on each row j of the mode-k loading matrix estimator, with maximum order of tensor time series as 3
#' @param k Mode of loading matrix.
#' @param D Eigenvalue matrix of sample covariance matrix, with dimension rk by rk.
#' @param Q Estimated mode-k loading matrix, with dimension Ik by rk.
#' @param C Estimated common component series, written in an array with dimension K+1 and mode-1 as the time mode.
#' @param Y Observed time series with missingness allowed, written in an array with dimension K+1 and mode-1 as the time mode.
#' @param j Integer representing the row of mode-k loading matrix. Value should be integers from minimum 1 to maximum Ik.
#' @param beta Lag parameter of the HAC type. Default is 0.
#'
#' @return A matrix of dimension rk by rk
#' @export
#'
#' @examples
#' K = 3;
#' TT = 10;
#' d = c(20,20,20);
#' r = c(2,2,2);
#' re = c(2,2,2);
#' eta = list(c(0,0), c(0,0), c(0,0));
#' coef_f = c(0.7, 0.3, -0.4, 0.2, -0.1);
#' coef_fe = c(-0.7, -0.3, -0.4, 0.2, 0.1);
#' coef_e = c(0.8, 0.4, -0.4, 0.2, -0.1);
#' data_test = tensor_gen(K,TT,d,r,re,eta, coef_f, coef_fe, coef_e);
#' data_miss = miss_gen(data_test$X);
#' data_est = miss_factor_est(data_miss, r);
#' D = diag(x=(svd(data_est$covMatrix[[2]])$d)[1:2], nrow=2, ncol=2);
#' sigmaD(2, D, data_est$A[[2]], data_est$imputation, data_miss, 2, 2);
#'
#'
#'
sigmaD <- function(k, D, Q, C, Y, j, beta=0){
  if (length(dim(C)) != length(dim(Y))){
    message('The tensor orders of C and Y are inconsistent.')
    return()
  } else if ( prod(dim(C)) != prod(dim(Y)) ){
    message('The dimensions of C and Y are inconsistent.')
    return()
  }
  
  K <- length(dim(C)) - 1
  TT <- dim(C)[1]
  if ( (k < 1)|(k > K) ){
    message('The input mode k is invalid with the order K of the given series.')
    return()
  }
  if ( (beta < 0)|(beta > (TT-1)) ){
    message('Parameter beta is either negative or too large.')
    return()
  }
  
  
  if (K==1){
    Sigma_HAC <- K1_D_nu(k+1, D, Q, dim(C), c(C), c(Y), j-1, 0)
    if (beta > 0){
      for (nu in 1:beta){
        D_nu <- K1_D_nu(k+1, D, Q, dim(C), c(C), c(Y), j-1, nu)
        D_nu <- (1- (nu/(1+beta))) * (t(D_nu) + D_nu)
        Sigma_HAC <- Sigma_HAC + D_nu
      }
    }
    return(Sigma_HAC)
    
  } else if (K==2){
    Sigma_HAC <- K2_D_nu(k+1, D, Q, dim(C), c(C), c(Y), j-1, 0)
    if (beta > 0){
      for (nu in 1:beta){
        D_nu <- K2_D_nu(k+1, D, Q, dim(C), c(C), c(Y), j-1, nu)
        D_nu <- (1- (nu/(1+beta))) * (t(D_nu) + D_nu)
        Sigma_HAC <- Sigma_HAC + D_nu
      }
    }
    return(Sigma_HAC)
    
  } else if (K==3){
    Sigma_HAC <- K3_D_nu(k+1, D, Q, dim(C), c(C), c(Y), j-1, 0)
    if (beta > 0){
      for (nu in 1:beta){
        D_nu <- K3_D_nu(k+1, D, Q, dim(C), c(C), c(Y), j-1, nu)
        D_nu <- (1- (nu/(1+beta))) * (t(D_nu) + D_nu)
        Sigma_HAC <- Sigma_HAC + D_nu
      }
    }
    return(Sigma_HAC)
    
  } else {
    message('Order K larger than 3 is not supported.')
    return()
  }
}