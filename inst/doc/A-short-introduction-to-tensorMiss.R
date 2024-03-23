## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=F, message=FALSE, eval=TRUE-------------------------------------
library(tensorMiss)

## -----------------------------------------------------------------------------
example <- array(1:24, dim=c(3,4,2))

## -----------------------------------------------------------------------------
example[3,1,1] <- NA
print(example)

## -----------------------------------------------------------------------------
example.1 <- unfold(example, 1)
print(example.1)

## -----------------------------------------------------------------------------
refold(example.1, 1, dim(example))==example

## -----------------------------------------------------------------------------
ttm(example, matrix(1:4,nrow=2), 3)

## -----------------------------------------------------------------------------
K <- 3 #order 3 at each time
TT <- 20 #time length
d <- c(30,30,30) #spatial dimensions
r <- c(2,2,2) #rank of core tensor
re <- c(2,2,2) #rank of common component in error
eta <- list(c(0,0), c(0,0), c(0,0)) #strong factors
coef_f <- c(0.7, 0.3, -0.4, 0.2, -0.1) #AR(5) coefficients for core factor
coef_fe <- c(-0.7, -0.3, -0.4, 0.2, 0.1) #AR(5) coefficients for common component in error
coef_e <- c(0.8, 0.4, -0.4, 0.2, -0.1) #AR(5) coefficients for idiosyncratic part in error
data_test <- tensor_gen(K,TT,d,r,re,eta, coef_f, coef_fe, coef_e)

## -----------------------------------------------------------------------------
data_miss <- miss_gen(data_test$X)

## -----------------------------------------------------------------------------
est_result <- miss_factor_est(data_miss, r)
fle(est_result$A[[1]], data_test$A[[1]])

## -----------------------------------------------------------------------------
qMSE(c(data_test$C), c(est_result$imputation), length(c(data_test$C))) #rMSE
qMSE(c(data_test$C), c(est_result$imputation), 100) #q-MSE

## -----------------------------------------------------------------------------
# computing the covariance matrix estimator
r2 <- r[2]
A2 <- data_test$A[[2]]
beta <- floor(0.2*(TT^0.25 * (d[2])^0.25))
D2 <- diag(x=(svd(est_result$covMatrix[[2]])$d)[1:r2], nrow=r2, ncol=r2)
# HAC_cov: HAC-type covariance matrix estimator
HAC_cov <- sigmaD(2, D2, est_result$A[[2]], est_result$imputation, data_miss, 1, beta)

# computing the rotation matrix
Amk <- data_test$A[[3]] %x% data_test$A[[1]]
R_ast <- 0
for (t in 1:TT){
  R_ast <- R_ast + unfold(data_test$Ft[t,,,],2) %*% t(Amk) %*% Amk %*% t(unfold(data_test$Ft[t,,,],2))
}
R_ast <- A2 %*% R_ast %*% t(A2)
R_ast <- R_ast/TT
Z2 <- diag(x = diag(t(A2) %*% A2), nrow=r2, ncol=r2)
Q2 <- A2 %*% diag(x=diag(solve(Z2))^0.5, nrow=r2, ncol=r2)
# H2: rotation matrix
H2 <- solve(D2) %*% t(est_result$A[[2]]) %*% R_ast %*% Q2 %*% solve(t(Q2)%*% Q2)

## -----------------------------------------------------------------------------
HAC_cov.eigen <- eigen(HAC_cov)
HAC_cov.sqrt <- HAC_cov.eigen$vectors %*% diag(sqrt(HAC_cov.eigen$values)) %*% solve(HAC_cov.eigen$vectors)
A2_1 <- (solve(HAC_cov.sqrt) %*% D2) %*% (matrix(est_result$A[[2]], nrow=d[2], ncol=r2)[1,] - (H2 %*% Q2[1,]))
A2_1

