library(Matrix)

sample_inverse_wishart_from_X <- function(X, nu) {
  n <- nrow(X)
  p <- ncol(X)
  
  if (n < p) {
    warning("n < p. The resulting scale matrix is singular.")
  }
  
  df <- nu + p + 1
  
  # 1. Scale X
  Y <- sqrt(nu) * X
  
  # 2. Extract the upper-triangular R using "reduced" QR decomposition
  R <- qr.R(qr(Y))  # R is (p, p)
  
  # 3. Generate the Bartlett factor L_W for a standard Wishart W(I, df)
  L_W <- matrix(0, p, p)
  diagonals <- sqrt(rchisq(p, df = df - (0:(p - 1))))
  diag(L_W) <- diagonals
  
  tril_count <- p * (p - 1) / 2
  L_W[lower.tri(L_W)] <- rnorm(tril_count)
  
  # 4. Solve L_W @ Z = R (L_W is lower-triangular, R is upper-triangular)
  Z <- forwardsolve(L_W, R)
  
  # 5. The final Inverse-Wishart sample
  sample <- t(Z) %*% Z
  
  return(sample)
}

library(reticulate)

# Import Python modules (requires msprime and stdpopsim installed in Python)
stdpopsim <- import("stdpopsim")
msprime   <- import("msprime")
np        <- import("numpy")

sim_susie_love <- function(chrom = "chr22",
                           start = 20e6,
                           end   = 21e6,
                           nu    = 1000L,
                           n_gwas = NULL,
                           J     = 1000L,
                           h2    = 0.01,
                           seed  = 5L) {

  
  
  # --- 1. Simulate ancestry via msprime (called through reticulate) ----------
  species <- stdpopsim$get_species("HomSap")
  contig  <- species$get_contig(chrom, left = start, right = end)
  
  ts <- msprime$sim_ancestry(
    samples           = N,
    population_size   = 20000L,
    recombination_rate = contig$recombination_map,
    model = list(
      msprime$DiscreteTimeWrightFisher(duration = 100L),
      msprime$StandardCoalescent()
    )
  )
  
  # --- 2. Add mutations -------------------------------------------------------
  mts    <- msprime$sim_mutations(ts, rate = 1.29e-8, random_seed = 1L)
  G_hap  <- mts$genotype_matrix()
  
  # --- 3. Collapse haplotypes to diploid genotypes (N x variants) -------------
  #   Python: G = (G_hap[:, 0::2] + G_hap[:, 1::2]).T
  #   G_hap is variants x haplotypes in tskit; we want individuals x variants
  G <- t(G_hap[, seq(1, ncol(G_hap), by = 2)] +
           G_hap[, seq(2, ncol(G_hap), by = 2)])
  rm(G_hap)
  
  # --- 4. Filter to common variants (MAF 0.05–0.95) --------------------------
  maf         <- colMeans(G) / 2
  common_mask <- (maf > 0.05) & (maf < 0.95)
  n_common    <- sum(common_mask)
  
  if (n_common == 0L) {
    stop("No common variants found. Try a larger region.")
  }
  
  G <- G[, common_mask, drop = FALSE]
  
  # --- 5. Randomly select a contiguous window of J SNPs ----------------------
  set.seed(seed)
  rng <- function(a) sample.int(a, size = 1L) - 1L   # 0-indexed like Python
  
  J_all <- max(ncol(G), J)
  win_start <- sample.int(J_all - J, size = 1L)
  G <- G[, win_start:(win_start + J - 1L), drop = FALSE]
  
  # --- 6. Generate summary statistics -----------------------------------------
  X <- scale(G, center = TRUE, scale = FALSE)   # mean-centre columns
  N <- nrow(X)
  J <- ncol(X)
  
  R <- cov(G)
  R[is.nan(R)] <- 0
  
  R0 <- sample_inverse_wishart_from_X(X / sqrt(N - 1), nu)
  
  # --- 7. Choose two nearly-uncorrelated causal SNPs --------------------------
  num_causal      <- 2L
  first_causal_SNP  <- sample.int(J, size = 1L)
  second_causal_SNP <- which.min(R[first_causal_SNP, ]^2)
  causal_SNPs_loc   <- c(first_causal_SNP, second_causal_SNP)
  
  # --- 8. True effect sizes (scaled to target h2) ----------------------------
  beta_true <- numeric(J)
  beta_true[causal_SNPs_loc] <- runif(num_causal) + 0.1
  
  Xb   <- X %*% beta_true
  sc   <- h2 / mean(Xb^2)
  beta_true <- beta_true * sqrt(sc)
  
  # --- 9. Phenotype and marginal association z-scores -------------------------
  y <- X %*% beta_true + sqrt(1 - h2) * rnorm(N)
  # y <- (y - mean(y)) / sd(y)
  
  v <- crossprod(X, y) / N          # J x 1 marginal statistics
  
  list(
    v               = as.numeric(v),
    R               = R,
    R0              = R0,
    beta_true       = beta_true,
    causal_SNPs_loc = causal_SNPs_loc
  )
}




# ### Plot R and R0
# # Print mean squared difference
# cat(mean((R - R0)^2), "\n")
# 
# # Plot side-by-side heatmaps
# par(mfrow = c(1, 2))
# 
# image(R[1:100, 1:100],
#       col  = hcl.colors(256, "viridis"),
#       axes = T,
#       main = "In-sample LD mat")
# 
# image(R0[1:100, 1:100],
#       col  = hcl.colors(256, "viridis"),
#       axes = T,
#       main = "Out-of-sample LD mat")
# 
# par(mfrow = c(1, 1))
# 
# B = tcrossprod(beta_true) 
# 
# J <- dim(R0)[1]
# nu_q <- nu
# nu_0 <- nu
# 
# c_target <- diag(R)
# 
# d0_inv_sqrt <- 1 / sqrt(diag(R0))
# LD0 <- diag(d0_inv_sqrt) %*% R0 %*% diag(d0_inv_sqrt)
# d_sqrt <- sqrt(diag(R))
# R_tilde <- diag(d_sqrt) %*% LD0 %*% diag(d_sqrt)
# 
# R_bar = .99 * R_tilde + 0.01 * diag(c_target)
# 
# cat("Target diagonal c:\n")
# print(round(c_target, 4))
# 
# cat("\nRunning isolated R_q optimizer...\n")
# results <- optimize_Rq_fixed(B, R_bar, v, c_target, N, nu_q, nu_0)
# 
# cat("\n=== VERIFICATION ===\n")
# actual_diag <- diag(results$R_q)
# max_error <- max(abs(actual_diag - c_target))
# 
# cat(sprintf("L-BFGS-B Converged: %s (0 means success)\n", results$convergence == 0))
# cat(sprintf("Number of Eigendecompositions evaluated: %d\n", results$evaluations))
# cat(sprintf("Max absolute error on R_q diagonal constraint: %.2e\n", max_error))
# 
# if (max_error < 1e-4) {
#   cat("SUCCESS: The solver perfectly locked onto the constraint.\n")
# } else {
#   cat("WARNING: Constraint not met.\n")
# }
# 
# cat("\nVisual Proof (First 5 diagonal elements):\n")
# print(data.frame(
#   Target_c = c_target[1:5],
#   Actual_Diag = actual_diag[1:5],
#   Diff = (c_target - actual_diag)[1:5]
# ))
# 
# results$lambda
# results$Y[1:10]
# results$L[1:10]

