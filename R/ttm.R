
#' Mode k product with matrix
#'
#' @description Performing k-mode matrix product of a tensor to a matrix
#' @param ten A multi-dimensional array with the mode-k dimension m.
#' @param A A matrix with dimension n by m.
#' @param k An integer specifying the tensor mode to perform k-mode matrix product.
#'
#' @return A multi-dimensional array with the k mode dimension n
#' @export
#'
#' @examples
#' ttm(array(1:24,c(3,4,2)), matrix(1:4,nrow =2), 3);
#'
#'
#'
ttm <- function(ten, A, k){
  result <- A %*% unfold(ten, k)
  dim_res <- dim(ten)
  dim_res[k] <- nrow(A)
  result <- refold(result, k, dim_res)
  result
}