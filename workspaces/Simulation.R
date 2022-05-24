rm(list = ls(all.names = TRUE))
#library(NPrior)
# source("Nprior_eval.R")
# library(Rcpp)
# sourceCpp("Neuro_Linear.cpp")
# sourceCpp("Neuro_ReLU_Linear.cpp")
# sourceCpp("Neuro_Linear_custom.cpp")

# Example in NPrior package
# library(MASS)
# data(Boston)
# str(Boston)
# attach(Boston)
# X = cbind(crim,indus,nox,rm,age,dis,tax,ptratio,black,lstat)
# X = scale(X)
# y = medv
# y = y-mean(y)
# n = nrow(X)
# p = ncol(X)
# theta <- NPrior_run(X, y)$THETA

######################### Simulation Studies ##################

## Data Generation
library(mvnfast)

# Function of generating data, low dim case
Gendata_low <- function(nsce,n,p){
  data <- list()
  for (m in 1:nsce){
    # Generate covariate matrix X
    # X_i \sim N(0,\Sigma) (each row i of X)
    Sigma <- matrix(NA,nrow=p,ncol=p)
    for (l in 1:p){
      for (k in 1:p){
        Sigma[l,k] <- 0.7^(abs(l-k))
      }
    }
    X <- matrix(NA,nrow=n,ncol=p)
    for (i in 1:n){
      X[i,] <- rmvn(1,rep(0,p),Sigma)
    }
    # Generate coefficient vector theta
    theta <- rep(NA,(p/10))
    for (i in 1:(p/10)){
      theta[i] <- (as.numeric(runif(1)<0.5)*2-1)*0.2
    }
    theta <- c(theta,rep(0,(p-p/10)))
    theta <- matrix(theta,ncol=1)
    # Generate epsilon 
    eps <- matrix(NA,ncol=1,nrow=n)
    for (i in 1:n){
      eps[i] <- rnorm(1)
    }
    # Calculate y
    y <- X%*%theta + eps
    data[[m]] <- list(X=X,theta=theta,y=y)
    # X (n \times p) matrix, theta (p \times 1) vector, y (p \times 1) vector
    cat(m,"\t")
  }
  return(data)
}

# Function of generating data, high dim case
Gendata_high <- function(nsce,n,p){
  theta_values <- c(0.4,0.45,0.5,0.55,0.6)
  data <- list()
  for (m in 1:nsce){
    # Generate covariate matrix X
    # X_i \sim N(0,\Sigma) (each row i of X)
    Sigma <- matrix(NA,nrow=p,ncol=p)
    for (l in 1:p){
      for (k in 1:p){
        Sigma[l,k] <- 0.7^(abs(l-k))
      }
    }
    X <- matrix(NA,nrow=n,ncol=p)
    for (i in 1:n){
      X[i,] <- rmvn(1,rep(0,p),Sigma)
    }
    # Generate coefficient vector theta
    theta <- rep(NA,5)
    for (i in 1:5){
      theta[i] <- (as.numeric(runif(1)<0.5)*2-1)*theta_values[i]
    }
    theta <- c(theta,rep(0,(p-5)))
    theta <- matrix(theta,ncol=1)
    # Generate epsilon 
    eps <- matrix(NA,ncol=1,nrow=n)
    for (i in 1:n){
      eps[i] <- rnorm(1)
    }
    # Calculate y
    y <- X%*%theta + eps
    data[[m]] <- list(X=X,theta=theta,y=y)
    # X (n \times p) matrix, theta (p \times 1) vector, y (p \times 1) vector
    cat(m,"\t")
  }
  return(data)
}

# Set seed
set.seed(4)
# number of scenarios
nsce <- 100

# data_low1 <- Gendata_low(nsce=nsce,n=100,p=50)
data_low2 <- Gendata_low(nsce=nsce,n=400,p=100)
data_low3 <- Gendata_low(nsce=nsce,n=200,p=50)
# data_high1 <- Gendata_high(nsce=nsce,n=100,p=300)
# data_high2 <- Gendata_high(nsce=nsce,n=150,p=1000)

# save(data_low1, file = "data_low1.RData")
save(data_low2, file = "data_low2.RData")
save(data_low3, file = "data_low3.RData")
# save(data_high1, file = "data_high1.RData")
# save(data_high2, file = "data_high2.RData")

# ##################### Simulation ################################
#
# # library(Matrix)
#
# # Function of calculating the cosine of the angle between true theta and thetahat (Cos)
# Cosfunc <- function(theta,thetahat){
#   Cos <- sum(theta*thetahat)/sqrt(sum(theta^2)*sum(thetahat^2))
#   return(Cos)
# }
#
# # Simulation
# # N-SpSL-L(Exact)
# # low dim 1
# set.seed(1234)
# M <- 100 # number of replications
# n <- 100
# p <- 50
# y <- data_low1[[1]]$y
# theta <- data_low1[[1]]$theta
# X <- data_low1[[1]]$X
# MSE <- rep(NA,M)
# MSEsd <- rep(NA,M)
# Cos <- rep(NA,M)
# for (m in 1:M){
#   theta_MCMC <- NPrior_run(X, y)$THETA
#   nMCMC <- ncol(theta_MCMC)
#   MSE_m <- rep(NA,nMCMC)
#   for (j in 1:nMCMC){
#     MSE_m[j] <- sum((theta-theta_MCMC[,j])^2)
#   }
#   MSE[m] <- mean(MSE_m)
#   MSEsd[m] <- sd(MSE_m)
#   thetahat <- apply(theta_MCMC,1,mean)
#   Cos[m] <- Cosfunc(theta,thetahat)
# }
# cat("MSE = ", mean(MSE),"\n MSE sd = ",mean(MSEsd),"\n Cos = ", mean(Cos))
#
#
#
#
#
