
#' Tensor refolding
#'
#' @description Performing to matrices tensorisation, which is the inverse process of unfolding
#' @param unfolding A matrix.
#' @param k An integer specifying the mode of array to refold from.
#' @param dim_vec A vector specifying the expected dimension of output array.
#'
#' @return A multi-dimensional array
#' @export
#' @import rTensor
#'
#' @examples
#' refold(matrix(1:9,nrow=3), 1, c(3,1,3));
#'
#'
refold <- function(unfolding, k, dim_vec){
  if (requireNamespace('rTensor', quietly = TRUE)){
    #a <- rTensor::as.tensor(unfolding)
    #a <- rTensor::k_fold(a, k, dim_vec)
    a <- as.tensor(unfolding)
    a <- k_fold(a, k, dim_vec)
    return(a@data)
  }
}