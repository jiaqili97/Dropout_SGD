####### random sample (y,X) ########
data_gen_linear <- function(n,d){
  set.seed(779)
  X <- matrix(rnorm(n*d, mean = 0, sd = 1), ncol=d)
  beta <- seq(0, 1, length.out = d + 1)
  beta <- beta[2:(d+1)]
  # beta <- c(1,2,5)
  epsilon <- rnorm(n, mean = 0, sd = 1)
  y <- X %*% beta + epsilon
  dat <- data.frame(y=y, x=X)
  return(list(dat,beta))
}

####### dropout matrix D ###########
dropout_mat <- function(n,d,p){
  set.seed(461)
  D_diag <- rbinom(n*d,1,p)
  D_list <- list()
  for (k in 1:n) {
    D_list[[k]] <- diag(D_diag[((k-1)*d+1):(k*d)]) 
  }
  return(D_list)
}

####### stochastic gradient with dropout ########
loss_grad_linear <- function(x, y, beta, D){
  loss_grad <- 2*D%*%x%*%(t(x) %*% D %*% beta - y) # transpose x
}


n <- 200000 # number of observations
d <- 10 # number of predictors
p <- 0.9 # probability of retaining the coordinate

data_linear <- data.matrix(data_gen_linear(n,d)[[1]])
X <- data_linear[, -1]
y <- data_linear[, 1]
beta <- data_gen_linear(n,d)[[2]]
D_list <- dropout_mat(n,d,p)


alpha_step <- 0.01

beta_sgd <- matrix(0, nrow = d, ncol = n)
beta_sgd_bar <- matrix(0, nrow = d, ncol = n)

beta_sgd[,1] <- 0 # can be changed to other initial points
beta_sgd_bar[,1] <- 0

## ASGD dropout estimation
for (k in 2:n) {
  beta_sgd[,k] <- beta_sgd[,k-1] - alpha_step*loss_grad_linear(X[k-1,], y[k-1], beta_sgd[,k-1], D_list[[k-1]])
  beta_sgd_bar[,k] <- ((k-1)*beta_sgd_bar[,k-1] + beta_sgd[,k])/k
}

# Plot convergence trace for each coordinate of beta
matplot(1:n, t(beta_sgd_bar), type = 'l', col = 1:d, 
        xlab = "n steps", 
        ylab = "estimated beta", 
        main = "Convergence of ASGD Dropout Iterates")
legend("topright", legend = paste("Beta", 1:d), col = 1:d, lty = 1)



# non-overlapping blocks {B_m}
M <- round(sqrt(n))
eta <- rep(0,M)
for (m in 1:M){
  eta[m] <- m^2
}
B <- list()
for (m in 1:(M-1)){
  B[[m]] <- seq.int(eta[m], eta[m+1]-1, by = 1)
}


################### Online estimation of long-run covariance #############################


psi <- numeric(n)
R_n <- matrix(0, nrow = d, ncol = n)
delta <- numeric(n)
K_n <- numeric(n)
H_n <- matrix(0, nrow = d, ncol = n)
Vcur_n <- array(0,dim=c(d,d,n))
V_n <- array(0,dim=c(d,d,n))
Sigma2_hat <- array(0,dim=c(d,d,n))


psi[1] <- 1

R_n[,1] <- beta_sgd[,1]
delta[1] <- 1
K_n[1] <- 1
H_n[,1] <- 1
Vcur_n[,,1] <- beta_sgd[,1] %o% beta_sgd[,1]
V_n[,,1] <- Vcur_n[,,1] + K_n[1]*beta_sgd_bar[,1]%o%beta_sgd_bar[,1] - 2*H_n[,1]%o%beta_sgd_bar[,1]
Sigma2_hat[,,1] <- V_n[,,1]

###### only needed for CI length ######
CI_length <- numeric(n)
Sigma2_hat_diag <- array(0, dim = c(d,d,n))
Sigma2_hat[,,1] <- diag(Sigma2_hat[,,1])
unit_v <- rep(1,d)/sqrt(d) # unit-length vector

for (k in 2:n){
  phi = (psi[k-1] + 1)^2
  
  if (k < phi){
    R_n[,k] <- R_n[,k-1] + beta_sgd[,k]
    psi[k] <- psi[k-1] 
    delta[k] <- delta[k-1] + 1
    K_n[k] <- K_n[k-1] - delta[k-1]^2 + delta[k]^2
    H_n[,k] <- H_n[,k-1] - delta[k-1]*R_n[,k-1] + delta[k]*R_n[,k]
    Vcur_n[,,k] <- Vcur_n[,,k-1] - R_n[,k-1]%o%R_n[,k-1] + R_n[,k]%o%R_n[,k]
  } else {
    R_n[,k] <- beta_sgd[,k]
    delta[k] <- 1
    psi[k] <- psi[k-1] + 1
    K_n[k] <- K_n[k-1] + 1
    H_n[,k] <- H_n[,k-1] + R_n[,k]
    Vcur_n[,,k] <- Vcur_n[,,k-1] + sum(R_n[,k]%o%R_n[,k])
  }
  
  V_n[,,k] <- Vcur_n[,,k] + K_n[k]*beta_sgd_bar[,k]%o%beta_sgd_bar[,k] - 2*H_n[,k]%o%beta_sgd_bar[,k]
  Sigma2_hat[,,k] <- V_n[,,k]/k
  
  ###### only needed for CI length with large dimension d ######
  Sigma2_hat_diag[,,k] <- diag(diag(Sigma2_hat[,,k]))
  
}
inner_sqrt_diag <- apply(Sigma2_hat_diag, 3, function(x) unit_v%*%x%*%unit_v)
CI_length <- 2*1.96*sqrt(abs(inner_sqrt_diag)/n)


# Extract diagonals of long-run covariance matrices, i.e., long-run variances
sigma2_hat <- matrix(0, nrow = d, ncol = n) 
for (j in 1:d) {
  sigma2_hat[j,] <- Sigma2_hat[j,j,]
}


# Plot convergence trace for long-run variances of coordinate of beta
matplot(1:n, t(sigma2_hat), type = 'l', col = 1:d, 
        xlab = "n steps", 
        ylab = "estimated long-run var", 
        main = "Estimated long-run variances of ASGD dropout")
# legend("topright", legend = paste("Variance", 1:d), col = 1:d, lty = 1)




# Plot convergence trace for long-run covariances of coordinate of beta
matplot(1:n, t(Sigma2_hat[1,,]), type = 'l', col = 1:d, 
        xlab = "n steps", 
        ylab = "estimated long-run cov", 
        main = "Estimated long-run covariances of ASGD dropout")
# legend("topright", legend = paste("Covariance", 1:d), col = 1:d, lty = 1)



# Plot CI length of 1d projection
plot(1:n,CI_length, type = 'l', ylab = 'Length of CI',
     xlab = 'n steps', main = 'Length of 95% CI', col = 'black')






