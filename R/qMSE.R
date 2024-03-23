
#' Quantile relative squared error
#'
#' @description Computing the q-quantile relative squared error as a generalised error measure on relative mean squared error
#' @param x_true True values, written in a vector of length n.
#' @param x_est Imputed or estimated values, written in a vector of length n.
#' @param q Number of partition intervals. If q equals n, then output is essentially relative mean squared error. Default is 100.
#'
#' @return A numeric number
#' @export
#'
#' @examples
#' qMSE(c(2, 3, 7, 1), c(-2, 0.5, 8, 2), 1);
#'
#'
#'
qMSE <- function(x_true, x_est, q=100){
  if (length(x_true) != length(x_est)){
    message('Lengths of the two vectors are different.')
    return()
  }
  if (q <= 0){
    message('Parameter q needs to be positive.')
    return()
  }
  return(partition_MSE(x_true, x_est, q))
}