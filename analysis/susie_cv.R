# The point of this script is to illustrate a simple idea for choosing
# L, the number of "single effect regressions", via a
# cross-validation-like approach.
library(matrixStats)
library(susieR)
library(ggplot2)
library(cowplot)
set.seed(1)

# First, simulate a data set using real genotype data.
X <- readRDS("../data/Thyroid.FMO2.1Mb.RDS")$X
storage.mode(X) <- "double"

# Take a random subset of 4,000 SNPs.
m <- ncol(X)
j <- sort(sample(m,4000))
X <- X[,j]

# Remove SNPs that show no variation.
j <- which(colSds(X) > 0)
X <- X[,j]

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
p  <- ncol(X)
j  <- sample(j,p1)

# Simulate b.
b        <- rep(0,p)
names(b) <- colnames(X)
b[j]     <- sample(c(-1,1),p1,replace = TRUE)

# Now simulate y.
e <- rnorm(n,sd = 1)
y <- drop(X %*% b + e)
y <- y - mean(y)
y <- y/sd(y)
y_train <- y[i0]
y_test  <- y[i1]
train   <- list(X = X_train,y = y_train)
test    <- list(X = X_test,y = y_test)
rm(X,X_train,X_test)
rm(y,y_train,y_test)

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
fit_susie <- susie(train$X,train$y,L = 10,standardize = TRUE,
                   min_abs_corr = 0,estimate_prior_method = "EM",
                   prior_tol = 0,check_null_threshold = -1,
                   verbose = TRUE)

# Run susie_rss on the training data.
n <- nrow(train$X)
fit_susie_rss <- susie_rss(zhat,R,n,L = 10,min_abs_corr = 0,
                           estimate_prior_method = "EM",
                           prior_tol = 0,check_null_threshold = -1,
                           verbose = TRUE)

# Run susie_rss on the training data, with out-of-sample LD (R_out).
fit_susie_rss_Rout <- susie_rss(zhat,R_out,n,L = 10,min_abs_corr = 0,
                                estimate_prior_method = "EM",prior_tol = 0,
                                check_null_threshold = -1,verbose = TRUE)

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

stop()

# Run susie on the training set with different settings of L, then
# predict the unseen Ys.
mse <- rep(0,10)
mse_rss <- rep(0,10)
for (L in 1:10) {
  input_ss <- compute_suff_stat(train$X,train$y,standardize = TRUE)
  ss <- univariate_regression(train$X,train$y)
  R <- with(input_ss,cov2cor(XtX))
  zhat <- with(ss,betahat/sebetahat)
  input_ss_test <- compute_suff_stat(test$X,test$y,standardize = TRUE)
  ss_test <- univariate_regression(test$X,test$y)
  R_test <- with(input_ss_test,cov2cor(XtX))
  zhat_test <- with(ss_test,betahat/sebetahat)
  fit <-  susie(train$X,train$y,L = L,standardize = TRUE,
                min_abs_corr = 0,estimate_prior_method = "EM",
                prior_tol = 0,check_null_threshold = -1,verbose = TRUE)
  fit_rss    <- susie_rss(zhat,R,n,min_abs_corr = 0,
                          estimate_prior_method = "EM",prior_tol = 0,
                          check_null_threshold = -1,verbose = TRUE)
  ypred      <- predict(fit,test$X)
  R_test[is.nan(R_test)] <- 0
  ypred_rss  <- predict(fit_rss,sqrt(n)*R_test)
  mse[L]     <- mean((test$y - ypred)^2)
  i <- which(!is.na(zhat_test))
  mse_rss[L] <- mean((zhat_test[i] - ypred_rss[i])^2)
}

par(mfcol = c(1,2))
plot(1:10,mse,type = "l",lwd = 1.5,col = "darkblue")
plot(1:10,mse_rss,type = "l",lwd = 1,col = "tomato",lty = "dashed")
