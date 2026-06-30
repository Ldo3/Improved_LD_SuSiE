# The point of this script is to illustrate the use of CV for model
# selection in susie, e.g., for choosing L.
library(matrixStats)
library(susieR)
library(ggplot2)
library(cowplot)
set.seed(1)

# First, simulate a data set using real genotype data.
X <- readRDS("../data/Thyroid.FMO2.1Mb.RDS")$X
storage.mode(X) <- "double"

# Remove SNPs with MAF < 1%
maf <- colMeans(X)/2
j <- which(maf > 0.01)
X <- X[,j]

# Take a random subset of 4,000 SNPs (just so the example does not
# take too long to run).
p <- ncol(X)
j <- sort(sample(p,4000))
X <- X[,j]
X <- scale(X,center = TRUE,scale = FALSE)

# Split the data into training and test sets, and remove any columns
# of X that show no variation in either the test or training sets.
n       <- nrow(X)
n_train <- 500
n_test  <- n - n_train
i0      <- seq(1,n_train)
i1      <- n_train + seq(1,n_test)
X_train <- X[i0,]
X_test  <- X[i1,]
j       <- which(colSds(X_train) > 0 & colSds(X_test) > 0)
X       <- X[,j]
X_train <- X_train[,j]
X_test  <- X_test[,j]

# Select the 3 causal SNPs.
p1 <- 3
p <- ncol(X)
j <- sample(1:p,p1)

# Simulate b.
b <- rep(0,p)
names(b) <- colnames(X)
b[j] <- sample(c(-1,1),p1,replace = TRUE)

# Now simulate y.
n <- nrow(X)
e <- rnorm(n,sd = 0.3) 
y <- drop(X %*% b + e)
y <- y - mean(y)
y_train <- y[i0]
y_test  <- y[i1]
train   <- list(X = X_train,y = y_train)
test    <- list(X = X_test,y = y_test)

stop()

# Compute the "in-sample" LD matrix using the training data, and
# compute summary statistics (association test z-scores).
out  <- with(train,compute_suff_stat(X,y,standardize = TRUE))
R    <- cov2cor(out$XtX)
out  <- with(train,univariate_regression(X,y))
zhat <- with(out,betahat/sebetahat)

# Compute an "out-of-sample" LD matrix using the test set. (Ideally
# we should use different data, but this will suffice for the purposes
# of this illustration.)
out   <- with(test,compute_suff_stat(X,y,standardize = TRUE))
R_out <- cov2cor(out$XtX)

# Run susie on the training data.
fit_susie <- susie(train$X,train$y,L = 5,standardize = TRUE,
                   min_abs_corr = 0,estimate_prior_method = "EM",
                   prior_tol = 0,check_null_threshold = -1,
                   verbose = TRUE)

# Run susie_rss on the training data.
n <- nrow(train$X)
fit_susie_rss <- susie_rss(zhat,R,n,L = 5,min_abs_corr = 0,
                           estimate_prior_method = "EM",
                           prior_tol = 0,check_null_threshold = -1,
                           verbose = TRUE)

# Run susie_rss on the training data, with out-of-sample LD (R_out).
fit_susie_rss_Rout <- susie_rss(zhat,R_out,n,L = 5,min_abs_corr = 0,
                                estimate_prior_method = "EM",prior_tol = 0,
                                check_null_threshold = -1,verbose = TRUE)

# Compare the coefficient estimates.
b_susie <- with(fit_susie,colSums(alpha * mu)/X_column_scale_factors)
b_susie_rss <- with(fit_susie_rss,colSums(alpha * mu)/X_column_scale_factors)
plot(b,b_susie,pch = 20)
plot(b_susie,b_susie_rss,pch = 20)

stop()

# Compare the PIP plots.
pip_plot <- function (fit, causal_snps = NULL, max_cs_size = 100, 
                      title = "") {
  colors <- c("black","greenyellow","red","magenta","darkorange",
              "olivedrab","cyan","dodgerblue","gold","forestgreen",
              "darkred")
  n    <- length(fit$pip)
  pdat <- data.frame(pos = 1:n,pip = fit$pip,cs = 0,causal = FALSE)
  num_cs <- length(fit$sets$cs)
  for (i in 1:num_cs) {
    j <- fit$sets$cs[[i]]
    if (length(j) <= max_cs_size)
      pdat[j,"cs"] <- i
  }    
  pdat[causal_snps,"causal"] <- TRUE
  pdat <- transform(pdat,cs = factor(cs))
  return(ggplot(pdat,aes(x = pos,y = pip,color = cs)) +
         geom_point() +
         geom_point(data = subset(pdat,causal),
                    mapping = aes(x = pos,y = pip),
                    shape = 1,size = 3,inherit.aes = FALSE) +
         scale_color_manual(values = colors) +
         labs(x = "SNP",y = "PIP",color = "CS",title = title) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}
p1 <- pip_plot(fit_susie,causal_snps = which(b != 0),title = "susie")
p2 <- pip_plot(fit_susie_rss,causal_snps = which(b != 0),title = "susie_rss")
p3 <- pip_plot(fit_susie_rss_Rout,causal_snps = which(b != 0),
               title = "susie_rss, out-of-sample LD")
print(plot_grid(p1,p2,p3,nrow = 3,ncol = 1))

# Implement a form of cross-validation by running susie on the
# training set with different settings of L, then predict the
# regression outcomes in the test set.
maxL <- 5
mse <- matrix(0,3,maxL)
rownames(mse) <- c("susie","susie_rss","susie_rss_Rout")
for (L in seq(1,maxL)) {
  cat(L,"")
  fit_susie <- susie(train$X,train$y,L = L,standardize = TRUE,
                     min_abs_corr = 0,estimate_prior_method = "EM",
                     prior_tol = 0,check_null_threshold = -1,
                     verbose = FALSE)
  ypred <- predict(fit_susie,test$X)
  mse["susie",L] <- mean((test$y - ypred)^2)
  n <- nrow(train$X)
  fit_susie_rss <- susie_rss(zhat,R,n,L = L,min_abs_corr = 0,
                             estimate_prior_method = "EM",
                             prior_tol = 0,check_null_threshold = -1,
                             verbose = FALSE)
  ypred <- predict(fit_susie_rss,test$X)
  mse["susie_rss",L] <- mean((test$y - ypred)^2)
  fit_susie_rss_Rout <- susie_rss(zhat,R_out,n,L = L,min_abs_corr = 0,
                                  estimate_prior_method = "EM",prior_tol = 0,
                                  check_null_threshold = -1,verbose = FALSE)
  ypred <- predict(fit_susie_rss_Rout,test$X)
  mse["susie_rss_Rout",L] <- mean((test$y - ypred)^2)
}
cat("\n")

plot(1:maxL,mse["susie",],type = "l",lwd = 2,col = "darkblue",
     xlab = "L",ylab = "MSE")
plot(1:maxL,mse["susie_rss",],lwd = 2,col = "dodgerblue")
plot(1:maxL,mse["susie_rss_Rout",],lwd = 2,col = "darkorange")
