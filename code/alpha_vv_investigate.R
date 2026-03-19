library(susieR)
library(Matrix)


setwd("~/Documents/Improved_LD_SuSiE")
# setwd("~/SuSiE-covariance/Improved_LD_SuSiE")

source("code/R_algorithms.R")
gtex = readRDS("data/Thyroid_ENSG00000132855.rds")

num_causal_SNP = 1
num_reps = 200


# Remove SNPs with MAF < 0.01
maf = apply(gtex, 2, function(x) sum(x)/2/length(x))
X0 = gtex[, maf > 0.01]
dim(X0)
X = na.omit(X0)
dim(X)
snp_total = ncol(X0)
all_true_effects = rep(0, num_reps)
choice_alphas = rep(0, num_reps)
opt_alpha = rep(0, num_reps)
sa_eigenvals_insample = rep(0, num_reps)
sa_eigenvals_outsample = rep(0, num_reps)

for (seed in 1:num_reps){
  print(seed)
  set.seed(seed)
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

  true_causal_SNPs = sample(c(1:J), size = num_causal_SNP, replace = FALSE)
  true_effect = rnorm(num_causal_SNP)
  all_true_effects[seed] <- true_effect / (abs(true_effect) + 1)
  # all_true_effects[seed] <- true_effect
  beta[true_causal_SNPs] <- true_effect
  # plot(beta, pch=16, ylab='effect size')
  y <- X_in %*% beta + rnorm(n)
  y = scale(y)

  v = (t(X_in) %*% y) / (n-1)
  choice_alphas[seed] = (1 / max(abs(v)))^2

  alpha_list = seq(0.9, 5, length.out = 41)
  for (i in 1:length(alpha_list)){
    alpha = alpha_list[i]
    R_minus = R - alpha * tcrossprod(v)
    ret = eigen(R_minus)
    sa_eigen = min(ret$values)
    if (sa_eigen < -1e-6){
      opt_alpha[seed] = alpha - .1
      break
    }
  }
  if (opt_alpha[seed]==0) opt_alpha[seed] = 5


  R_minus = R - tcrossprod(v)
  ret = eigen(R_minus)
  sa_eigenvals_insample[seed] = min(ret$values)

  R_minus = R_hat - tcrossprod(v)
  ret = eigen(R_minus)
  sa_eigenvals_outsample[seed] = min(ret$values)
}

library(ggplot2)

# Create a data frame
df <- data.frame(x = all_true_effects, y = choice_alphas)

# ggplot scatter plot
ggplot(df, aes(x = x, y = y)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot",
       x = "true effect beta",
       y = "our choice of alpha = (1 / max(v))^2")


df <- data.frame(x = all_true_effects, y = opt_alpha)
# ggplot scatter plot
ggplot(df, aes(x = x, y = y)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot",
       x = "true effect beta",
       y = "optimal alpha")

df <- data.frame(x = all_true_effects, y = sa_eigenvals_insample)
# ggplot scatter plot
ggplot(df, aes(x = x, y = y)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot",
       x = "true effect beta",
       y = "smallest eigenvals of R_in - vvT")

df <- data.frame(x = all_true_effects, y = sa_eigenvals_outsample)
# ggplot scatter plot
ggplot(df, aes(x = x, y = y)) +
  geom_point(color = "blue") +
  labs(title = "Scatter Plot",
       x = "true effect beta",
       y = "smallest eigenvals of R_out - vvT")




