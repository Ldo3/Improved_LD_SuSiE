source("./code/R_algorithms.R")
gtex = readRDS("./data/Thyroid_ENSG00000132855.RDS")
maf = apply(gtex, 2, function(x) sum(x)/2/length(x))
X0 = gtex[, maf > 0.01]
X0 = na.omit(X0)
snp_total = ncol(X0)

p = 1000
num_causal_SNPs = 1
# Start from a random point on the genome
indx_start = sample(1: (snp_total - p), 1)
X = X0[, indx_start:(indx_start + p -1)]
# View(cor(X)[1:10, 1:10])

n = nrow(X)
## sub-sample into two
out_sample = sample(1:n, 100)
X_out = X[out_sample, ]
X_in = X[setdiff(1:n, out_sample), ]
sum(is.na(X_out))

rm_p = c(which(diag(cov(X_in))==0), which(diag(cov(X_out))==0))
length(rm_p)
indx_p = setdiff(1:p, rm_p)
X_in = X_in[, indx_p]
X_out = X_out[, indx_p]

## Standardize both sample matrices
X_in <- scale(X_in)
X_out <- scale(X_out)

## out-sample LD matrix
R_hat = cor(X_out)
R = cor(X_in)

## generate data from in-sample X matrix
p = ncol(X_in)
J = ncol(X_in)
beta <- rep(0,p)
n = nrow(X_in)
true_causal_SNP = sample(c(1:p), num_causal_SNPs, replace = F)
beta[true_causal_SNP] <- rnorm(num_causal_SNPs) * 2
# plot(beta, pch=16, ylab='effect size')
y <- X_in %*% beta + rnorm(n)
y = scale(y)

v = t(X_in) %*% y / (n-1)

errors_full = rep(0, 3)
times_full = rep(0, 3)
errors_RSpec = rep(0, 4)
times_RSpec = rep(0, 4)

## proj B
start_time <- proc.time()[3]
R1 <- proj_B(R_hat, v)
time_sec <- proc.time()[3] - start_time
times_full[1] <- time_sec
errors_full[1] <- sum((R1 - R)^2)

start_time <- proc.time()[3]
R1 <- project_psd_shifted_iterative(R_hat, v)$P
time_sec <- proc.time()[3] - start_time
times_RSpec[1] <- time_sec
errors_RSpec[1] <- sum((R1 - R)^2)

## proj B, diag 1
start_time <- proc.time()[3]
R2 <- proj_Dykstra(R_hat, v, max_iter = 200L, tol = 1e-4)$R
time_sec <- proc.time()[3] - start_time
times_full[2] <- time_sec
errors_full[2] <- sum((R2 - R)^2)

start_time <- proc.time()[3]
R2 <- dykstra_rank1_iterative_opt2(R_hat, list(v),
                                   tol = 1e-4,
                                   max_iter = 200L,
                                   k = 5L)$R
time_sec <- proc.time()[3] - start_time
times_RSpec[2] <- time_sec
errors_RSpec[2] <- sum((R2 - R)^2)

## multi-trait
v_list = list(v)
for (t in 1:10){
  beta_ = rep(0, J)
  casual_SNPs = sample(c(1:p), num_causal_SNPs, replace = F)
  true_effect = rnorm(num_causal_SNPs) * 2
  beta_[casual_SNPs] <- true_effect
  # plot(beta, pch=16, ylab='effect size')
  y_ <- X_in %*% beta_ + rnorm(n)
  y_ = scale(y_)
  v_ = (t(X_in) %*% y_) / (n-1)
  v_list <- c(v_list, list(v_))
}

start_time <- proc.time()[3]
R3 <- proj_Dykstra_mtraits(R_hat, v_list, max_iter = 200L, tol = 1e-4)$R
time_sec <- proc.time()[3] - start_time
times_full[3] <- time_sec
errors_full[3] <- sum((R3 - R)^2)

start_time <- proc.time()[3]
# ret <- dykstra_rank1_iterative_opt1(R_hat, v_list,
#                                     max_iter = 200L,
#                                     tol = 1e-4,
#                                     tol_eig = 1e-2,
#                                     k_initial = 5L) ## So slow bruh

ret <- dykstra_rank1_iterative_opt2(R_hat, v_list,
                                    tol = 1e-4,
                                    max_iter = 200L,
                                    k = 10L)
time_sec <- proc.time()[3] - start_time
R3 = ret$R
times_RSpec[3] <- time_sec
errors_RSpec[3] <- sum((R3 - R)^2)

# start_time <- proc.time()[3]
# R3 <- dykstra_rank1_iterative_opt2(R_hat, v_lists,
#                                    max_iter = 200L,
#                                    tol = 1e-4,
#                                    tol_eig = 1e-4)$R
# time_sec <- proc.time()[3] - start_time
# times_RSpec[4] <- time_sec
# errors_RSpec[4] <- sum((R3 - R)^2)

times_full
times_RSpec
errors_full
errors_RSpec

