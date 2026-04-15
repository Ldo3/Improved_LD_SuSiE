library(matrixStats)
library(susieR)
X <- readRDS("../data/Thyroid.FMO2.1Mb.RDS")$X
storage.mode(X) <- "double"

# Remove SNPs that show no variation.
j <- which(colSds(X) > 0)
X <- X[,j]

# Select the causal SNPs.
p  <- ncol(X)
p1 <- 3
j  <- sample(j,p1)

# Simulate b.
b        <- rep(0,p)
names(b) <- colnames(X)
b[j]     <- sample(c(-1,1),p1,replace = TRUE)

# Simulate y.
n <- nrow(X)
e <- rnorm(n,sd = 0.3)
y <- drop(X %*% b + e)
y <- y/sd(y)
causal_snps <- j

# Split into training and test sets.
n_train <- 500
n_test  <- n - n_train
i       <- seq(1,n_train)
j       <- n_train + seq(1,n_test)
train   <- list(X = X[i,],y = y[i])
test    <- list(X = X[j,],y = y[j])

# Run susie on the training set with different settings of L, then
# predict the unseen Ys.
mse <- rep(0,10)
for (L in 1:10) {
  fit <-  susie(train$X,train$y,L = L,standardize = FALSE,
                min_abs_corr = 0,estimate_prior_method = "EM",
                prior_tol = 0,verbose = TRUE)
  ypred <- predict(fit,test$X)
  mse[L] <- mean((test$y - ypred)^2)
}
