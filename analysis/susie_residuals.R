# remotes::install_github("stephenslab/susieR@output_residuals")
library(susieR)
library(matrixStats)
print(packageVersion("susieR"))
set.seed(1)
p1      <- 2
min_maf <- 0.03
n_out   <- 100

# Load the genotype data.
X <- readRDS("../data/Thyroid.FMO2.1Mb.RDS")$X
storage.mode(X) <- "double"
X <- X[,1:1000]

# Filter the SNPs by MAF.
maf <- colMeans(X)/2
maf <- pmin(maf,1 - maf)
j   <- which(maf > min_maf)
X   <- X[,j]
print(dim(X))

# Choose the causal SNPs.
p <- ncol(X)
j <- sample(p,p1)

# Simulate b.
b        <- rep(0,p)
names(b) <- colnames(X)
b[j]     <- sample(c(-1,1),p1,replace = TRUE)

# Simulate y.
n <- nrow(X)
e <- rnorm(n,sd = 0.85)
y <- drop(X %*% b + e)
y <- y/sd(y)

# Compute the "in-sample" LD matrix and an "out-of-sample" LD matrix.
R <- cor(X)
i <- sample(n,n_out)
R_out <- cor(X[i,])

# Compute the summary statistics.
out  <- univariate_regression(X,y)
zhat <- with(out,betahat/sebetahat)

# Run susie_rss with the "in-sample" LD matrix.
fit <- susie_rss(zhat,R,L = 4,n = n,estimate_prior_method = "EM",
                 tol = 1e-6,min_abs_corr = 0,verbose = TRUE)
par(mfrow = c(1,2))
susie_plot(fit,y = "PIP",b = (b != 0))

# Run susie_rss with the "out-of-sample" LD matrix.
fit_out <- susie_rss(zhat,R_out,L = 4,n = n,estimate_prior_method = "EM",
                     tol = 1e-6,min_abs_corr = 0,verbose = TRUE)
susie_plot(fit_out,y = "PIP",b = (b != 0))

fit_out_L2 <- susie_rss(zhat,R_out,L = 2,n = n,estimate_prior_method = "EM",
                        tol = 1e-6,min_abs_corr = 0,verbose = TRUE)

susie_plot(fit,y = "PIP",b = (b != 0))
susie_plot(fit_out_L2,y = "PIP",b = (b != 0))

hist(fit$XtR,breaks = 64)
hist(fit_out_L2$XtR,breaks = 64)
