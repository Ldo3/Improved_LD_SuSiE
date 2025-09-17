library(susieR)
library(Matrix)    
library(RSpectra)

setwd("~/Documents/Improved_LD_SuSiE")
source("code/R_algorithms.R")
gtex = readRDS("data/Thyroid_ENSG00000132855.rds")

seed = 16
set.seed(seed)
print(paste("This is seed number", seed))

# Remove SNPs with MAF < 0.01
maf = apply(gtex, 2, function(x) sum(x)/2/length(x))
X0 = gtex[, maf > 0.01]
dim(X0)
X = na.omit(X0)
dim(X)
snp_total = ncol(X0)
n = nrow(X0)
J = 200
# Start from a random point on the genome
indx_start = sample(1: (snp_total - J), 1)
X = X0[, indx_start:(indx_start + J -1)]
# View(cor(X)[1:10, 1:10])

## sub-sample into two
out_sample = sample(1:n, 100)
X_out = X[out_sample, ]
X_in = X[setdiff(1:n, out_sample), ]
sum(is.na(X_out))

rm_p = c(which(diag(cov(X_in))==0), which(diag(cov(X_out))==0))
length(rm_p)
indx_p = setdiff(1:J, rm_p)
X_in = X_in[, indx_p]
X_out = X_out[, indx_p]

## Standardize both sample matrices
X_in <- scale(X_in)
X_out <- scale(X_out)

## out-sample LD matrix
R_hat = cor(X_out)
R = cor(X_in)
# View(R_hat[1:10, 1:10])
# View(R[1:10, 1:10])
# View(cor(X)[1:10, 1:10])

## generate data from in-sample X matrix
J = ncol(X_in)
beta <- rep(0, J)
n = nrow(X_in)
casual_SNPs = c(50, 100, 150)
true_effect = c(2, -3, 2)
beta[casual_SNPs] <- true_effect
# plot(beta, pch=16, ylab='effect size')
y <- X_in %*% beta + rnorm(n)
y = scale(y)

## compute summary statistics 
sumstats <- univariate_regression(X_in, y)
z_scores <- sumstats$betahat / sumstats$sebetahat
# susie_plot(z_scores, y = "z", b=beta)
min_cor = 0.5
# L_fit = 3
L_fit = 10
## 1. fit the susie-rss model with in-sample R
fitted_rss1 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = R, var_y = var(y), L = L_fit,
                         estimate_residual_variance = TRUE,
                         min_abs_corr=min_cor)
# summary(fitted_rss1)$cs
# p1 = susie_plot(fitted_rss1, y="PIP", b=beta)

## 2. fit the model with out-sample R
fitted_rss2 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = R_hat, var_y = var(y), L = L_fit,
                         estimate_residual_variance = FALSE,
                         min_abs_corr=min_cor)
# will have problem non-positive cov if estimate_residual_variance = TRUE
# summary(fitted_rss2)$cs
# p2 = susie_plot(fitted_rss2, y="PIP", b=beta) ## miss the true or does not run


## 3. adjusted by identity matrix
lambda = 0.1
R_hat_lambd = (1-lambda) * R_hat + lambda * diag(J)
fitted_rss3 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = R_hat_lambd, var_y = var(y), L = L_fit,
                         estimate_residual_variance = F,
                         min_abs_corr=min_cor)
# will have problem non-positive cov if estimate_residual_variance = TRUE
# summary(fitted_rss3)$cs
# susie_plot(fitted_rss3, y="PIP", b=beta) 

## 4. using our method + Dykstra projection
v = (t(X_in) %*% y) / (n-1)

ret = proj_Dykstra(R=R_hat, v=v)
R_hat_plus_diag1 = ret$R
info = ret$info
fitted_rss4 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = R_hat_plus_diag1, var_y = var(y), L = L_fit,
                         estimate_residual_variance = F,
                         min_abs_corr=min_cor)


## 5. Only project on the space \{R: R\geq vv^T\}
R_hat_minus = R_hat - tcrossprod(v)
eigen_R = eigen(R_hat_minus)
eigen_R$values

V <- eigen_R$vectors
D_plus <- diag(pmax(eigen_R$values, 0))

R_hat_plus <- V %*% D_plus %*% t(V) + tcrossprod(v)

fitted_rss5 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = R_hat_plus, var_y = var(y), L = L_fit,
                         estimate_residual_variance = F,
                         min_abs_corr=min_cor)
# summary(fitted_rss5)$cs
# susie_plot(fitted_rss5, y="PIP", b=beta)

## 6. vv^T -> set diag = 1
v = (t(X_in) %*% y) / (n-1)

ret = proj_Dykstra(R=matrix(0, nrow=J, ncol=J), v=v)
lst_vvT_diag1 = ret$R
info0 = ret$info
fitted_rss6 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = lst_vvT_diag1, var_y = var(y), L = L_fit,
                         estimate_residual_variance = F,
                         min_abs_corr=min_cor)

## 7. R = vv^T

fitted_rss7 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = tcrossprod(v), var_y = var(y), L = L_fit,
                         estimate_residual_variance = F,
                         min_abs_corr=min_cor)

## 8. Identity matrix
fitted_rss8 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                         R = diag(1, J), var_y = var(y), L = L_fit,
                         estimate_residual_variance = F,
                         min_abs_corr=min_cor)

par(mfrow = c(3, 3))
susie_plot(z_scores, y = "z", b=beta)
title('Z-score')
susie_plot(fitted_rss1, y="PIP", b=beta)
title('in-sample LD')
susie_plot(fitted_rss7, y="PIP", b=beta)
title('R = vv^T')

susie_plot(fitted_rss2, y="PIP", b=beta)
title('out-sample LD')
susie_plot(fitted_rss3, y="PIP", b=beta)
title('Regularized lambda=0.1')
susie_plot(fitted_rss8, y="PIP", b=beta)
title('Identity matrix (lambda=1)')

susie_plot(fitted_rss4, y="PIP", b=beta)
title('Dykstra projection(R_hat)')
susie_plot(fitted_rss5, y="PIP", b=beta)
title('Proj_{vv}(R_hat)')
susie_plot(fitted_rss6, y="PIP", b=beta)
title('Dykstra projection(0)')


