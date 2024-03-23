
#' Tensor unfolding
#'
#' @description Performing to multi-dimensional arrays tensor unfolding, also known as matricization
#' @param ten A multi-dimensional array.
#' @param k An integer specifying the mode of array to unfold.
#'
#' @return A matrix
#' @export
#' @import rTensor
#'
#' @examples
#' unfold(array(1:24, dim=c(3,4,2)), 2);
#'
#'
#'
unfold <- function(ten,k){
  if (requireNamespace('rTensor', quietly = TRUE)){
    #a <- rTensor::as.tensor(ten)
    #a <- rTensor::rs_unfold(a, k)
    a <- as.tensor(ten)
    a <- rs_unfold(a, k)
    return(a@data)
  }
}