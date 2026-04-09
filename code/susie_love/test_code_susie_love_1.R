setwd("~/Documents/Improved_LD_SuSiE/code/susie_love")
source("data-generate-coalescent.R")
source("susie-rss-my.R")
source("susie_love_ver_fixed_lambda.R")

# ---- 1. Simulate data from coalescent
seed  <- 1

T_split = 25
chrom <- "chr22"
start <- 20e6
end   <- 21e6
N     <- 20000L
N0    <- 200L
J     <- 1000L
h2    <- 0.01

res <- sim_2pop(chrom, start, end, 
                T_split = T_split, 
                n_gwas = N, 
                n_ref = N0, 
                J = J, 
                h2 = h2, 
                seed = seed) 

v               <- res$v
R               <- res$R
R0              <- res$R0
beta_true       <- res$beta_true
causal_SNPs_loc <- res$causal_SNPs_loc

true_indices = causal_SNPs_loc

plot(v, pch = 16)
points(true_indices, v[true_indices], pch = 16, col = "red")


# ---- 2. Fit susie with in-sample covariance matrix
library(susieR)
fit <- susie_suff_stat(
  XtX = N * R, 
  Xty = N * v, 
  n = N,
  yty = N
)

susie_plot(fit, y = "PIP", b = beta_true, main = "In-sample cov mat")

# ---- 3. Fit susie with out-of-sample covariance matrix
fit_R0 <- susie_suff_stat(
  XtX = N * R0, 
  Xty = N * v, 
  n = N,
  yty = N, 
  L = 5
)

susie_plot(fit_R0, y = "PIP", b = beta_true, main = "Out-of-sample cov mat")


## adjust the diagonal and fit again
d0_inv_sqrt <- 1 / sqrt(diag(R0))
LD0 <- diag(d0_inv_sqrt) %*% R0 %*% diag(d0_inv_sqrt)
d_sqrt <- sqrt(diag(R))
R_tilde <- diag(d_sqrt) %*% LD0 %*% diag(d_sqrt)
fit_R_tilde <- susie_suff_stat(
  XtX = N * R_tilde, 
  Xty = N * v, 
  n = N,
  yty = N
)
susie_plot(fit_R_tilde, y = "PIP", b = beta_true, main = "Out-of-sample cov mat (corrected diag)")


fit_R_tilde_L5 <- susie_suff_stat(
  XtX = N * R_tilde, 
  Xty = N * v, 
  n = N,
  yty = N, 
  L = 5
)
susie_plot(fit_R_tilde_L5, y = "PIP", b = beta_true, main = "Out-of-sample cov mat (corrected diag), L = 5")



# ---- 4. Fit susie-love
L = 5
k <- 200  ## rank of Z
nu_q <- 1000
nu_0 <- 1000

J <- dim(R0)[1]
c_diag <- diag(R)

lambd = 0.001
R_bar = (1 - lambd) * R_tilde + lambd * diag(c_diag) 

# PRECOMPUTE: log_det_Rbar never changes, so do it once outside the loop
log_det_Rbar <- as.numeric(determinant(R_bar, logarithm = TRUE)$modulus)

## Init Z for R_q 
z_current <- initialize_Z_from_prior_R(R_bar, nu_0, nu_q, c_diag, k)

## Init susie parameters
alpha <- matrix(1 / J, nrow = L, ncol = J)
mu <- matrix(0, nrow = L, ncol = J)
sigma2_post = matrix(0.01, nrow = L, ncol = J)

b_res = list(alpha = alpha, 
             mu = mu,
             sigma2_post = sigma2_post, 
             sigma02 = rep(0.1, L))

max_iters = 30
for (iter in 1:max_iters) {
  # 1. --- Update b ---  
  Z <- matrix(z_current, nrow = J, ncol = k)
  s <- 1 + rowSums(Z^2)
  alpha_diag <- sqrt(c_diag / s)
  F_mat <- Z * alpha_diag
  D_diag <- alpha_diag^2
  R_q = F_mat %*% t(F_mat) + diag(D_diag)
  
  b_res = susie_rss_my(b_res, v, R_q, N, max_iter = 20, tol = 1e-4)
  
  # Update B
  alpha = b_res$alpha
  mu = b_res$mu
  sigma2_post = b_res$sigma2_post
  sigma02 = b_res$sigma02
  B_current = compute_E_bbT(alpha, mu, sigma2_post)
  
  # 2. --- Update R_q ---
  scalar_V <- (nu_q + J + 1) / nu_q
  V_matrix <- scalar_V * (nu_0 * R_bar + N * tcrossprod(v))
  a_const <- (nu_0 + J + 2) / (2 * J)
  b_const <- N / (2 * J)
  c_const <- 1 / (2 * J)
  
  optim_res <- optim(
    par = z_current,
    fn = objective_Rq,
    gr = gradient_Rq,
    B = B_current,          
    V = V_matrix,           
    c_diag = c_diag, J = J, k = k,
    a_const = a_const, b_const = b_const, c_const = c_const,
    method = "L-BFGS-B",
    control = list(maxit = 30) 
  )
  
  z_current <- optim_res$par
  
  Z_updated <- matrix(z_current, nrow = J, ncol = k)
  
  
  # 3. --- Update nu_q and nu_0 ---
  # Reconstruct fast Rq and Rq_inv components using the newly updated Z
  
  s_updated <- 1 + rowSums(Z_updated^2)
  alpha_updated <- sqrt(c_diag / s_updated)
  
  F_mat_updated <- Z_updated * alpha_updated
  D_diag_updated <- alpha_updated^2
  D_inv_diag_updated <- 1 / D_diag_updated
  
  F_T_Dinv_F <- t(F_mat_updated) %*% (F_mat_updated * D_inv_diag_updated)
  M <- diag(k) + F_T_Dinv_F
  
  chol_M <- chol(M)
  log_det_M <- 2 * sum(log(diag(chol_M)))
  log_det_Rq <- sum(log(D_diag_updated)) + log_det_M
  
  M_inv <- chol2inv(chol_M)
  D_inv_F <- F_mat_updated * D_inv_diag_updated
  Rq_inv <- diag(D_inv_diag_updated) - D_inv_F %*% tcrossprod(M_inv, D_inv_F)
  
  # Precompute the traces for the nu objective function
  trace_Rq_inv_Rbar <- sum(Rq_inv * R_bar)
  trace_Rq_inv_vvT <- as.numeric(crossprod(v, Rq_inv %*% v))
  
  nu_optim_res <- optim(
    par = c(nu_q, nu_0),
    fn = objective_nu,
    J = J, N = N,
    log_det_Rbar = log_det_Rbar,
    log_det_Rq = log_det_Rq,
    trace_Rq_inv_Rbar = trace_Rq_inv_Rbar,
    trace_Rq_inv_vvT = trace_Rq_inv_vvT,
    method = "L-BFGS-B",
    lower = c(1, 1),  # Bounded to keep the Inv-Wishart well-defined
    control = list(maxit = 5)
  )
  
  # Extract updated parameters
  nu_q <- nu_optim_res$par[1]
  nu_0 <- nu_optim_res$par[2]
  
  
  # 4. --- Print Status ---
  current_loss <- compute_loss(b_res, Z_updated, nu_q, nu_0, v, N, R_bar, c_diag)
  cat(sprintf("Iter %d | Loss: %.4f | nu_q: %.1f | nu_0: %.1f\n", 
              iter, current_loss, nu_q, nu_0))
}

Z <- matrix(z_current, nrow = J, ncol = k)
s <- 1 + rowSums(Z^2)
alpha_diag <- sqrt(c_diag / s)
F_mat <- Z * alpha_diag
D_diag <- alpha_diag^2
R_q = F_mat %*% t(F_mat) + diag(D_diag)

cs_results <- get_cs(alpha = b_res$alpha, R = R_q, coverage = 0.95)
# Generate the official-looking plot

susie_plot_alpha(
  alpha = b_res$alpha, 
  R = cov2cor(R_q),               
  b = beta_true, 
  main = "SuSiE_LovE plot"
)


print(paste0("Distance from learned cov mat is ", sum((R_q - R)^2)))
print(paste0("Distance from out-of-sample cov mat is ", sum((R0 - R)^2)))
print(paste0("Distance from adjusted diag cov mat is ", sum((R_tilde - R)^2)))


