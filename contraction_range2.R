N <- 500

n <- 100
d <- 5


X <- matrix(rnorm(n * d), n, d)

eigen_max_X <- max(eigen(t(X)%*%X)$values)

alpha_max <- 2/eigen_max_X
alpha <- 0.0154


p <- 0.9
A_sq <- list()
for (s in 1:N) {
  D <- diag(rbinom(d,1,p), nrow = d)
  I <- diag(1,nrow = d)
  A <- I - alpha*D%*%t(X)%*%X%*%D
  A_sq[[s]] <- t(A)%*%A
}
A_sq_mean <- Reduce("+", A_sq) / length(A_sq)



eigen_max_A_sq <- max(eigen(A_sq_mean)$values)


alpha_max
alpha
eigen_max_A_sq < 1




