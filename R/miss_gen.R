
#' Assignment of missingness to tensor time series
#'
#' @description Assign missingness to a given order-K tensor time series, where the maximum K is 4
#' @param dt Tensor time series, written in an array with dimension K+1 and mode-1 as the time mode.
#' @param type Type of missingness, where "random" is random missing with probability p, "simul" is missingness on the last half along all dimensions, "mix" is a mixture of "random" and "simul". Default is "random".
#' @param p If type is "random", then each entry is randomly missing with probability 1-p. Default is 0.7.
#'
#' @return A multi-dimensional array with dimension K+1, where mode-1 is the time mode and missing entries are denoted by NA
#' 
#' @export
#' @import stats
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
#' data_test = tensor_gen(K,TT,d,r,re,eta, coef_f, coef_fe, coef_e);
#' miss_gen(data_test$X);
#'
#'
miss_gen <- function(dt, type="random", p=0.7){
  if (requireNamespace('stats', quietly = TRUE)){
    
    if (length(dim(dt)) > 5){
      message('Order of the given tensor is too large.')
      return()
    } else if (length(dim(dt)) == 1){
      message('Scalar time series is not supported.')
      return()
    }
    
    if ((p<0)|(p>1)){
      message('Probability p is invalid.')
      return()
    }
    
    
    
    if (type == "random"){
      miss_array <- array(rbinom(prod(dim(dt)), 1, p), dim=dim(dt))
      miss_array[which(miss_array==0)] <- NA
      return(dt * miss_array)
      
    } else if ( (type == "simul")|(type == "mix") ) {
      if (type == "simul"){
        miss_array <- array(rbinom(prod(dim(dt)), 1, 1), dim=dim(dt))
      } else {
        miss_array <- array(rbinom(prod(dim(dt)), 1, p), dim=dim(dt))
      }
      
      if (length(dim(dt)) == 2){
        TT <- dim(dt)[1]
        I1 <- dim(dt)[2]
        miss_array[( (floor(0.5*TT)) + 1 ):TT, ( (floor(0.5*I1)) + 1 ):I1] <- 0
        miss_array[which(miss_array==0)] <- NA
        return(dt * miss_array)
        
      } else if (length(dim(dt)) == 3){
        TT <- dim(dt)[1]
        I1 <- dim(dt)[2]
        I2 <- dim(dt)[3]
        miss_array[( (floor(0.5*TT)) + 1 ):TT, ( (floor(0.5*I1)) + 1 ):I1, ( (floor(0.5*I2)) + 1 ):I2] <- 0
        miss_array[which(miss_array==0)] <- NA
        return(dt * miss_array)
        
      } else if (length(dim(dt)) == 4){
        TT <- dim(dt)[1]
        I1 <- dim(dt)[2]
        I2 <- dim(dt)[3]
        I3 <- dim(dt)[4]
        miss_array[( (floor(0.5*TT)) + 1 ):TT, ( (floor(0.5*I1)) + 1 ):I1, ( (floor(0.5*I2)) + 1 ):I2, ( (floor(0.5*I3)) + 1 ):I3] <- 0
        miss_array[which(miss_array==0)] <- NA
        return(dt * miss_array)
        
      } else if (length(dim(dt)) == 5){
        TT <- dim(dt)[1]
        I1 <- dim(dt)[2]
        I2 <- dim(dt)[3]
        I3 <- dim(dt)[4]
        I4 <- dim(dt)[5]
        miss_array[( (floor(0.5*TT)) + 1 ):TT, ( (floor(0.5*I1)) + 1 ):I1, ( (floor(0.5*I2)) + 1 ):I2, ( (floor(0.5*I3)) + 1 ):I3, ( (floor(0.5*I4)) + 1 ):I4] <- 0
        miss_array[which(miss_array==0)] <- NA
        return(dt * miss_array)
        
      }
      
    } else {
      message('Type parameter is invalid.')
      return()
    } 
    
  }
}