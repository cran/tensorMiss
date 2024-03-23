
#' Factor loading error
#'
#' @description Computing the column space distance between two matrix
#' @param A1 A matrix of m rows and n columns.
#' @param A2 A matrix of m rows and l columns where l can equal n.
#'
#' @return A numeric number
#' @export
#'
#' @examples
#' fle(matrix(1:12, nrow=4), matrix(11:22, nrow=4));
#'
#'
#'
fle <- function(A1, A2){
  rot.A1 <- qr.Q(qr(A1))%*%t(qr.Q(qr(A1)))
  rot.A2 <- qr.Q(qr(A2))%*%t(qr.Q(qr(A2)))
  return( norm(rot.A1 - rot.A2, type = '2') )
}