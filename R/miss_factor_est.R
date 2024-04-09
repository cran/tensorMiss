
#' Estimation of tensor factor models with missing data
#'
#' @description Estimate the factor structure on an order-K tensor at each time t, with maximum K as 3 and missing entries allowed
#' @param dt Tensor time series, written in an array with dimension K+1 and mode-1 as the time mode.
#' @param r Rank of core tensors, written in a vector of length K. First value as 0 is to denote unknown rank which would be automatically estimated using ratio-based estimators. Default is 0.
#' @param delta Non-negative number as the correction parameter for rank estimation. Default is 0.2.
#'
#' @return A list containing the following:
#' r: a vector representing either the given rank or the estimated rank, with length K;
#' A: a list of estimated K factor loading matrices;
#' Ft: the estimated core factor series, as multi-dimensional array with dimension K+1, where mode-1 is the time mode;
#' imputation: the imputed common component time series, as multi-dimensional array with dimension K+1, where mode-1 is the time mode;
#' covMatrix: a list of estimated covariance matrix which are used to estimate loading matrices;
#' 
#' 
#' @export
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @useDynLib tensorMiss, .registration = TRUE
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
#' miss_factor_est(data_miss, r);
#'
#'
#'
miss_factor_est <- function(dt, r=0, delta=0.2){
  if (requireNamespace('RcppEigen', quietly = TRUE)){
    if (length(dim(dt)) > 4){
      message('Order of the given tensor is too large. If there is need to extend to higher order, check the source Rcpp file or consult the author.')
      return()
    } else if (length(dim(dt)) == 1){
      message('Scalar time series is not supported.')
      return()
    }
    
    K <- length(dim(dt)) - 1
    
    if ( (r[1]!=0)&(K!=length(r)) ){
      message('Length of r is inconsistent with the order of given tensor time series.')
      return()
    }
    if (delta < 0){
      message('Parameter delta needs to at least 0.')
      return()
    }
    
    # start estimation -----------------------------------------------------------
    if (K==1){
      TT <- dim(dt)[1]
      I1 <- dim(dt)[2]
      
      # loading matrix estimation
      est_loading <- list()
      for (k in 2:length(dim(dt)) ){
        est_loading[[k-1]] <- K1_cov_est(c(dt), dim(dt), k)
      }
      
      # rank estimation if necessary
      if ( (r[1])==0 ){
        xi_1 <- delta * I1*(TT^(-0.5) + I1^(-0.5))
        e1_val <- svd( est_loading[[1]] )$d[1:floor(I1/2)] + xi_1
        r1 <- which.min( ((e1_val[-1])/(e1_val[-length(e1_val)])) )
      } else {
        r1 <- r[1] 
      } # rank estimation end here
      
      A1.est <- matrix(svd( est_loading[[1]] )$u[, (1:r1)], nrow=I1, ncol=r1)
      
      # factor series estimation
      F.est <- array(0, dim=c(TT, r1))
      for (t in 1:TT ){
        F_t.est <- K1_Ft_est(c(dt[t,]), A1.est)
        F.est[t,] <- array(F_t.est, dim=c(r1))
      }
      
      # imputation
      Y.imp <- ttm(F.est,  A1.est, 2)
      return(list(r=c(r1), A=list(A1.est), Ft=F.est, imputation=Y.imp, covMatrix=est_loading))
      
    } else if (K==2){
      TT <- dim(dt)[1]
      I1 <- dim(dt)[2]
      I2 <- dim(dt)[3]
      
      # loading matrix estimation
      est_loading <- list()
      for (k in 2:length(dim(dt)) ){
        est_loading[[k-1]] <- K2_cov_est(c(dt), dim(dt), k)
      }
      
      # rank estimation if necessary
      if ( (r[1])==0 ){
        xi_1 <- delta * I1*I2*((TT*I2)^(-0.5) + I1^(-0.5))
        e1_val <- svd( est_loading[[1]] )$d[1:floor(I1/2)] + xi_1
        r1 <- which.min( ((e1_val[-1])/(e1_val[-length(e1_val)])) )
        
        xi_2 <- delta * I1*I2*((TT*I1)^(-0.5) + I2^(-0.5))
        e2_val <- svd( est_loading[[2]] )$d[1:floor(I2/2)] + xi_2
        r2 <- which.min( ((e2_val[-1])/(e2_val[-length(e2_val)])) )
      } else {
        r1 <- r[1] 
        r2 <- r[2]
      } # rank estimation end here
      
      A1.est <- matrix(svd( est_loading[[1]] )$u[, (1:r1)], nrow=I1, ncol=r1)
      A2.est <- matrix(svd( est_loading[[2]] )$u[, (1:r2)], nrow=I2, ncol=r2)
      
      # factor series estimation
      F.est <- array(0, dim=c(TT, r1, r2))
      for (t in 1:TT ){
        F_t.est <- K2_Ft_est(c(dt[t,,]), A1.est, A2.est)
        F.est[t,,] <- array(F_t.est, dim=c(r1, r2))
      }
      
      # imputation
      Y.imp <- ttm(F.est,  A1.est,2)
      Y.imp <- ttm(Y.imp, A2.est,3)
      return(list(r=c(r1,r2), A=list(A1.est,A2.est), Ft=F.est, imputation=Y.imp, covMatrix=est_loading))
      
    } else if (K==3){
      TT <- dim(dt)[1]
      I1 <- dim(dt)[2]
      I2 <- dim(dt)[3]
      I3 <- dim(dt)[4]
      
      # loading matrix estimation
      est_loading <- list()
      for (k in 2:length(dim(dt)) ){
        est_loading[[k-1]] <- K3_cov_est(c(dt), dim(dt), k)
      }
      
      # rank estimation if necessary
      if ( (r[1])==0 ){
        xi_1 <- delta * I1*I2*I3*((TT*I2*I3)^(-0.5) + I1^(-0.5))
        e1_val <- svd( est_loading[[1]] )$d[1:floor(I1/2)] + xi_1
        r1 <- which.min( ((e1_val[-1])/(e1_val[-length(e1_val)])) )
        
        xi_2 <- delta * I1*I2*I3*((TT*I1*I3)^(-0.5) + I2^(-0.5))
        e2_val <- svd( est_loading[[2]] )$d[1:floor(I2/2)] + xi_2
        r2 <- which.min( ((e2_val[-1])/(e2_val[-length(e2_val)])) )
        
        xi_3 <- delta * I1*I2*I3*((TT*I2*I1)^(-0.5) + I3^(-0.5))
        e3_val <- svd( est_loading[[3]] )$d[1:floor(I3/2)] + xi_3
        r3 <- which.min( ((e3_val[-1])/(e3_val[-length(e3_val)])) )
      } else {
        r1 <- r[1] 
        r2 <- r[2]
        r3 <- r[3]
      } # rank estimation end here
      
      A1.est <- matrix(svd( est_loading[[1]] )$u[, (1:r1)], nrow=I1, ncol=r1)
      A2.est <- matrix(svd( est_loading[[2]] )$u[, (1:r2)], nrow=I2, ncol=r2)
      A3.est <- matrix(svd( est_loading[[3]] )$u[, (1:r3)], nrow=I3, ncol=r3)
      
      # factor series estimation
      F.est <- array(0, dim=c(TT, r1, r2, r3))
      for (t in 1:TT ){
        F_t.est <- K3_Ft_est(c(dt[t,,,]), A1.est, A2.est, A3.est)
        F.est[t,,,] <- array(F_t.est, dim=c(r1, r2, r3))
      }
      
      # imputation
      Y.imp <- ttm(F.est,  A1.est,2)
      Y.imp <- ttm(Y.imp, A2.est,3)
      Y.imp <- ttm(Y.imp, A3.est,4)
      return(list(r=c(r1,r2,r3), A=list(A1.est,A2.est,A3.est), Ft=F.est, imputation=Y.imp, covMatrix=est_loading))
    }
  }
}