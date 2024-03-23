
#' Data generation of tensor time series with factor structure
#'
#' @description Generate an order-K tensor at each time t, with the first mode as the time mode and maximum allowed K is 4
#' @param K Order of the generated tensor at each time t.
#' @param TT Length of time series.
#' @param d Dimensions of each mode of the tensor, written in a vector of length K.
#' @param r Rank of the core tensors, written in a vector of length K.
#' @param re Rank of the cross-sectional common error core tensors, written in a vector of length K.
#' @param eta Quantities controlling factor strengths in each factor loading matrix, written in a list of K vectors.
#' @param coef_f AR(5) coefficients for the factor series, written in a vector of length 5.
#' @param coef_fe AR(5) coefficients for the common component in error series, written in a vector of length 5.
#' @param coef_e AR(5) coefficients for the idiosyncratic component in error series, written in a vector of length 5.
#' @param heavy_tailed Whether to generate data from heavy-tailed distribution. If FALSE, generate from N(0,1); if TRUE, generate from t-distribution. Default is FALSE.
#' @param t_df The degree of freedom for t-distribution if heavy_tailed = TRUE. Default is 3.
#' @param seed Random seed required for reproducibility. Default is 2023.
#'
#'
#' @return A list containing the following:
#' X: the generated tensor time series, as multi-dimensional array with dimension K+1, where mode-1 is the time mode;
#' A: a list of K factor loading matrices;
#' C: the generated common component time series, as multi-dimensional array with dimension K+1, where mode-1 is the time mode;
#' Ft: the generated core factor series, as multi-dimensional array with dimension K+1, where mode-1 is the time mode;
#' 
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
#' tensor_gen(K,TT,d,r,re,eta, coef_f, coef_fe, coef_e);
#'
#'
#'
tensor_gen <- function(K,TT,d,r,re,eta, coef_f, coef_fe, coef_e, heavy_tailed=FALSE, t_df=3, seed = 2023){
  if ((K<1)|(K>4)){
    message('K can only be integers in 1,2,3,4.')
    return()
  }
  
  Data_test = data_gen(K,TT,d,r,re,eta, coef_f, coef_fe, coef_e, heavy_tailed, t_df, seed)
  
  if (K==1){
    #construct common component
    F_series <- array(Data_test$F_ts, dim=c(TT, r))
    C_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_series[t,] <- array(Data_test$Ao %*% c(F_series[t,]), dim = d)
    }
    #construct error common component
    E_common <- array(Data_test$E_common, dim=c(TT, re))
    C_E_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_E_series[t,] <- array(Data_test$Aeo %*% c(E_common[t,]), dim = d)
    }
    #construct full data
    X_series <- array(Data_test$E_idio, dim=c(TT, d))
    X_series <- X_series + C_series + C_E_series
    return(list(X=X_series, A=Data_test$A, C=C_series, Ft=F_series))
    
  } else if (K==2){
    #construct common component
    F_series <- array(Data_test$F_ts, dim=c(TT, r))
    C_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_series[t,,] <- array(Data_test$Ao %*% c(F_series[t,,]), dim = d)
    }
    #construct error common component
    E_common <- array(Data_test$E_common, dim=c(TT, re))
    C_E_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_E_series[t,,] <- array(Data_test$Aeo %*% c(E_common[t,,]), dim = d)
    }
    #construct full data
    X_series <- array(Data_test$E_idio, dim=c(TT, d))
    X_series <- X_series + C_series + C_E_series
    return(list(X=X_series, A=Data_test$A, C=C_series, Ft=F_series))
    
  } else if (K==3){
    #construct common component
    F_series <- array(Data_test$F_ts, dim=c(TT, r))
    C_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_series[t,,,] <- array(Data_test$Ao %*% c(F_series[t,,,]), dim = d)
    }
    #construct error common component
    E_common <- array(Data_test$E_common, dim=c(TT, re))
    C_E_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_E_series[t,,,] <- array(Data_test$Aeo %*% c(E_common[t,,,]), dim = d)
    }
    #construct full data
    X_series <- array(Data_test$E_idio, dim=c(TT, d))
    X_series <- X_series + C_series + C_E_series
    return(list(X=X_series, A=Data_test$A, C=C_series, Ft=F_series))
    
  } else if (K==4){
    #construct common component
    F_series <- array(Data_test$F_ts, dim=c(TT, r))
    C_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_series[t,,,,] <- array(Data_test$Ao %*% c(F_series[t,,,,]), dim = d)
    }
    #construct error common component
    E_common <- array(Data_test$E_common, dim=c(TT, re))
    C_E_series <- array(0, dim=c(TT,d))
    for (t in 1:TT){
      C_E_series[t,,,,] <- array(Data_test$Aeo %*% c(E_common[t,,,,]), dim = d)
    }
    #construct full data
    X_series <- array(Data_test$E_idio, dim=c(TT, d))
    X_series <- X_series + C_series + C_E_series
    return(list(X=X_series, A=Data_test$A, C=C_series, Ft=F_series))
    
  }
}