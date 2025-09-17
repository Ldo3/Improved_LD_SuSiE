library(susieR)
setwd("~/Documents/Improved_LD_SuSiE")
source("code/SuSiE_rss.R")
source("code/R_algorithms.R")


# ## simulate multi casual SNP
# n = 500
# J = 200
# set.seed(1)
# X = matrix(rnorm(n*J), nrow = n, ncol = J)
# num_casual_SNP = 3
# gamma = sample(c(1:J), size = num_casual_SNP, replace = FALSE)
# b = rnorm(num_casual_SNP) * 3
# beta = rep(0, J)
# beta[gamma] = b
# y = X %*% beta + rnorm(n)
# print(gamma)
# ## standardize and compute sufficient statistics
# X = scale(X) / sqrt(n-1)
# y = scale(y) / sqrt(n-1)
# V_xx = t(X) %*% X
# V_xy = t(X) %*% y
# 
# ## start training model
# sigma2 = 1
# sigma02 = 1
# L = 3
# max_iter = 2
# 
# library(ggplot2)
# library(patchwork)
# 
# # plots <- list()
# 
# 
# b_bar = matrix(0, nrow = L, ncol = J)
# b_bar2 = matrix(0, nrow = L, ncol = J)
# alphas = matrix(0, nrow = L, ncol = J)
# mus = matrix(0, nrow = L, ncol = J)
# sigmas = matrix(0, nrow = L, ncol = J)
# 
# 
# for (iter in (1:max_iter)){
#   par(mfrow = c(3, 3))
#   V_xy_bar = V_xy - V_xx %*% colSums(b_bar)  ## residual signal
#   
#   for (ell in 1:L){
#     V_xy_bar_ell = V_xy_bar + V_xx %*% b_bar[ell, ]  ## add back ell-th signal
#     
#     susie_plot(V_xy_bar_ell, y = "z", b=beta)
#     
#     ret = SER(V_xx, V_xy_bar_ell, n, sigma2, sigma02)
#     alphas[ell, ] = ret$alpha
#     mus[ell, ] = ret$mus
#     sigmas[ell, ] = ret$sigma12
#     b_bar[ell, ] = alphas[ell, ] * mus[ell, ]
#     b_bar2[ell, ] = alphas[ell, ] * (mus[ell, ]^2 + sigmas[ell, ])
#     V_xy_bar = V_xy_bar_ell - V_xx %*% b_bar[ell, ]
#     
#     susie_plot(ret$alpha, y = "PIP", b=beta)
#     susie_plot(V_xy_bar, y = "z", b=beta)
#   }
# }


###############################################
## using gtex data
setwd("~/Documents/Improved_LD_SuSiE")
source("code/R_algorithms.R")
gtex = readRDS("data/Thyroid_ENSG00000132855.rds")

## seed present: 1
seed = 2
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


## generate data from in-sample X matrix
J = ncol(X_in)
beta <- rep(0, J)
n = nrow(X_in)
num_casual_SNP = 1
gamma = sample(c(1:J), size = num_casual_SNP, replace = FALSE)
b = rnorm(num_casual_SNP) * 3
beta = rep(0, J)
beta[gamma] = b
y = X_in %*% beta + rnorm(n)
y = scale(y) 

V_xy = t(X_in) %*% y / (n-1)

## 1. In-sample covariance matrix
V_xx = R
sigma2 = 1
sigma02 = 1
L = 3
max_iter = 10

b_bar = matrix(0, nrow = L, ncol = J)
b_bar2 = matrix(0, nrow = L, ncol = J)
alphas = matrix(0, nrow = L, ncol = J)
mus = matrix(0, nrow = L, ncol = J)
sigmas = matrix(0, nrow = L, ncol = J)


for (iter in (1:max_iter)){
  par(mfrow = c(3, 3))
  V_xy_bar = V_xy - V_xx %*% colSums(b_bar)  ## residual signal
  
  for (ell in 1:L){
    V_xy_bar_ell = V_xy_bar + V_xx %*% b_bar[ell, ]  ## add back ell-th signal
    
    susie_plot(V_xy_bar_ell, y = "z", b=beta)
    
    ret = SER(V_xx, V_xy_bar_ell, n, sigma2, sigma02)
    alphas[ell, ] = ret$alpha
    mus[ell, ] = ret$mus
    sigmas[ell, ] = ret$sigma12
    b_bar[ell, ] = alphas[ell, ] * mus[ell, ]
    b_bar2[ell, ] = alphas[ell, ] * (mus[ell, ]^2 + sigmas[ell, ])
    V_xy_bar = V_xy_bar_ell - V_xx %*% b_bar[ell, ]
    
    susie_plot(ret$alpha, y = "PIP", b=beta)
    susie_plot(V_xy_bar, y = "z", b=beta)
    
  }
  mtext(paste0("In-sample R, Iteration ", iter), outer = TRUE, cex = 1.5)
}

## 1. Out-sample covariance matrix
V_xx = R_hat
sigma2 = 1
sigma02 = 1
L = 3
max_iter = 10

b_bar = matrix(0, nrow = L, ncol = J)
b_bar2 = matrix(0, nrow = L, ncol = J)
alphas = matrix(0, nrow = L, ncol = J)
mus = matrix(0, nrow = L, ncol = J)
sigmas = matrix(0, nrow = L, ncol = J)


for (iter in (1:max_iter)){
  par(mfrow = c(3, 3))
  V_xy_bar = V_xy - V_xx %*% colSums(b_bar)  ## residual signal
  
  for (ell in 1:L){
    V_xy_bar_ell = V_xy_bar + V_xx %*% b_bar[ell, ]  ## add back ell-th signal
    
    susie_plot(V_xy_bar_ell, y = "z", b=beta)
    
    ret = SER(V_xx, V_xy_bar_ell, n, sigma2, sigma02)
    alphas[ell, ] = ret$alpha
    mus[ell, ] = ret$mus
    sigmas[ell, ] = ret$sigma12
    b_bar[ell, ] = alphas[ell, ] * mus[ell, ]
    b_bar2[ell, ] = alphas[ell, ] * (mus[ell, ]^2 + sigmas[ell, ])
    V_xy_bar = V_xy_bar_ell - V_xx %*% b_bar[ell, ]
    
    susie_plot(ret$alpha, y = "PIP", b=beta)
    susie_plot(V_xy_bar, y = "z", b=beta)
    
  }
  mtext(paste0("Out-sample R, Iteration ", iter), outer = TRUE, cex = 1.5)
}

## 3. Projected R
ret = proj_Dykstra(R=R_hat, v=V_xy)
V_xx = ret$R
sigma2 = 1
sigma02 = 1
L = 3
max_iter = 10

b_bar = matrix(0, nrow = L, ncol = J)
b_bar2 = matrix(0, nrow = L, ncol = J)
alphas = matrix(0, nrow = L, ncol = J)
mus = matrix(0, nrow = L, ncol = J)
sigmas = matrix(0, nrow = L, ncol = J)


for (iter in (1:max_iter)){
  par(mfrow = c(3, 3))
  V_xy_bar = V_xy - V_xx %*% colSums(b_bar)  ## residual signal
  
  for (ell in 1:L){
    V_xy_bar_ell = V_xy_bar + V_xx %*% b_bar[ell, ]  ## add back ell-th signal
    
    susie_plot(V_xy_bar_ell, y = "z", b=beta)
    
    ret = SER(V_xx, V_xy_bar_ell, n, sigma2, sigma02)
    alphas[ell, ] = ret$alpha
    mus[ell, ] = ret$mus
    sigmas[ell, ] = ret$sigma12
    b_bar[ell, ] = alphas[ell, ] * mus[ell, ]
    b_bar2[ell, ] = alphas[ell, ] * (mus[ell, ]^2 + sigmas[ell, ])
    V_xy_bar = V_xy_bar_ell - V_xx %*% b_bar[ell, ]
    
    susie_plot(ret$alpha, y = "PIP", b=beta)
    susie_plot(V_xy_bar, y = "z", b=beta)
    
  }
  mtext(paste0("Projected R, Iteration ", iter), outer = TRUE, cex = 1.5)
}


