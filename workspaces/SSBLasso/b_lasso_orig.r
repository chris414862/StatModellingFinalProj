# Package names

# library(mnormt) # rmnorm
# library(VGAM) # rinv.gaussian
# library(miscTools) # colMeans
packages <- c("mlbench", "mnormt", "VGAM", "miscTools", "monomvn", "SSLASSO")
## monomvn package - blasso function
## docs:  https://cran.r-project.org/web/packages/monomvn/monomvn.pdf

print(packages)


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
print("installed packages")
print(installed_packages)
if (any(installed_packages == FALSE)) {
  # install.packages(packages[!installed_packages], "~/R/x86_64-pc-linux-gnu-library/3.6")
  install.packages(packages[!installed_packages])#, "~/R/x86_64-pc-linux-gnu-library/3.6")
}
print("after")
print(installed_packages)

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


data(BostonHousing)
bh <- BostonHousing
# source("gibbsBLasso.r")
source("helper_funcs.r")
print("DATASET SUMMARY:")
summary(bh)
# set predictor variables
X <- bh[, !names(bh) %in% c("medv", "chas")]
X_dims <- dim(X)
n <- X_dims[1] # rows: Num data points
k <- X_dims[2] # cols: Num features

# set dependent variable
y <- bh$medv
# dim(y) <- c(n, 1)# give proper dims
# print("Starting 5 values from in X:")
# print(X[1:5,])
# print("Starting 5 values from in y:")
# print(y[1:5])

# Scale input
y = standard(y)
X = standard(X)
# print("Standardized values from in X:")
# print(X[1:5,])
# print(typeof(X))
# print(is.matrix(X)) # Should be true
# print("Standardized values from in y:")
# print(y[1:5])

lambda1 <- .1
length <- 50

lambda0 <- seq(lambda1, n,length=length)
post_data <- SSLASSO(X, y, lambda1=lambda1, lambda0=lambda0)#, RJ = FALSE)
betas = post_data$beta
## beta values for the MAP estimate of posterior when using diff values of lambda0

beta_dims = dim(betas)
beta_k = beta_dims[1]
beta_len = beta_dims[2]
stopifnot(k==beta_k, length==beta_len) # These values should be equal

reps = beta_dims[2]
for (rep in seq(1, reps, length=reps)){
    if (rep > 1){
        break
    }
    
    beta_vec = betas[,rep]
    dim(beta_vec) <- c(k,1) # give proper dims
    
    # print("AFTER")
    # print("Betas:")
    # print(beta_vec)
    #
    # print("Betas dims:")
    # print( dim(beta_vec))
    # print("Betas type:")
    # print( str(beta_vec))

}

print(paste("k:", k, "reps:", reps, sep=" "))


