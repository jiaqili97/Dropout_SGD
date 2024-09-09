N <- 500 # number of replications

n <- 100 # sample size
d <- 5 # dimension


X <- matrix(rnorm(n * d), n, d) # full design matrix

eigen_max_X <- max(eigen(t(X)%*%X)$values)

alpha_max <- 2/eigen_max_X # maximum allowed learning rate
alpha <- 0.0154 # you can change this input learning rate according to alpha_max


p <- 0.9 # success probability of dropout
A_sq <- list()
for (s in 1:N) {
  D <- diag(rbinom(d,1,p), nrow = d) # dropout matrix
  I <- diag(1,nrow = d)
  A <- I - alpha*D%*%t(X)%*%X%*%D # random coefficient matrix
  A_sq[[s]] <- t(A)%*%A
}
A_sq_mean <- Reduce("+", A_sq) / length(A_sq)



eigen_max_A_sq <- max(eigen(A_sq_mean)$values) # contraction contant


alpha_max # maximum allowed learning rate
alpha # your input learning rate
eigen_max_A_sq < 1 # if this line returns "TRUE", then the contraction holds




