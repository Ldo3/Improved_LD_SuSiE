# dykstra_rank1_iterative.R
# Dykstra projection onto intersection of shifted PSD cones C_t = { S : S - u_t u_t^T >= 0 }
# Uses RSpectra::eigs_sym to iteratively find negative eigenpairs and avoid full dense eigendecompositions.

library(Matrix)    # for sparse matrices if desired
library(RSpectra)  # for eigs_sym

symm <- function(A) {
  # enforce symmetry
  return((A + t(A)) / 2)
}

# Project matrix M onto set { S : S - u u^T >= 0 }
# using iterative eigensolver to find negative eigenpairs of A = M - u u^T.
project_psd_shifted_iterative <- function(M, u,
                                          max_k_fraction = 0.5,
                                          k_initial = 5L,
                                          k_step_factor = 2L,
                                          tol_eig = 1e-8,
                                          eigs_opts = list()) {
  # M : numeric matrix (n x n)
  # u : numeric vector length n
  n <- nrow(M)
  A <- symm(M - tcrossprod(u, u))  # A = M - u u^T
  max_k <- max(1L, as.integer(max_k_fraction * n))
  k <- min(max(1L, as.integer(k_initial)), n - 2L)
  used_dense <- FALSE
  
  # Try iterative loop
  repeat {
    # try RSpectra eigs_sym to get k smallest algebraic eigenvalues
    # RSpectra::eigs_sym uses which = "SA" for smallest algebraic by default
    success_iter <- TRUE
    try({
      ev <- eigs_sym(A, k = k, which = "SA", opts = eigs_opts)
    }, silent = TRUE) -> tc
    
    if (!exists("ev") || is.null(ev$values) || length(ev$values) != k) {
      success_iter <- FALSE
    }
    
    if (!success_iter) {
      # fallback to dense
      ev_full <- eigen(A, symmetric = TRUE)
      vals_f <- ev_full$values
      vecs_f <- ev_full$vectors
      vals_pos <- pmax(vals_f, 0)
      P_plus <- vecs_f %*% (vals_pos * t(vecs_f))  # vecs * diag(vals_pos) * t(vecs)
      P <- tcrossprod(u, u) + symm(P_plus)
      return(list(P = symm(P), used_dense = TRUE))
    }
    
    # sort eigenpairs ascending
    idx <- order(ev$values)
    vals <- ev$values[idx]
    vecs <- ev$vectors[, idx, drop = FALSE]
    
    # negative eigenvalues among returned
    neg_mask <- vals < -tol_eig
    neg_vals <- vals[neg_mask]
    if (length(neg_vals) == 0) {
      # likely no negative eigs -> A is PSD
      P_plus <- A
      P <- tcrossprod(u, u) + symm(P_plus)
      return(list(P = symm(P), used_dense = FALSE))
    }
    
    # If the largest of the returned eigenvalues is >= -tol_eig,
    # we've likely captured all negative eigenvalues (no smaller negatives remain).
    if (vals[length(vals)] >= -tol_eig) {
      # build correction C = sum_{lambda<0} lambda v v^T
      C <- matrix(0.0, n, n)
      neg_vecs <- vecs[, neg_mask, drop = FALSE]
      for (j in seq_len(ncol(neg_vecs))) {
        vj <- neg_vecs[, j, drop = FALSE]
        C <- C + as.numeric(neg_vals[j]) * (vj %*% t(vj))
      }
      P_plus <- A - C
      P <- tcrossprod(u, u) + symm(P_plus)
      return(list(P = symm(P), used_dense = FALSE))
    }
    
    # else, vals[k] < 0 meaning more negative eigenvalues may exist
    if (k >= max_k || k >= (n - 2L)) {
      # fallback to dense
      ev_full <- eigen(A, symmetric = TRUE)
      vals_f <- ev_full$values
      vecs_f <- ev_full$vectors
      vals_pos <- pmax(vals_f, 0)
      P_plus <- vecs_f %*% (vals_pos * t(vecs_f))
      P <- tcrossprod(u, u) + symm(P_plus)
      return(list(P = symm(P), used_dense = TRUE))
    }
    
    # increase k and retry
    k <- min(n - 2L, as.integer(k * k_step_factor))
    if (k <= 0L) k <- 1L
  } # end repeat
}

# Dykstra algorithm using the iterative projection above
dykstra_rank1_iterative <- function(R, u_list,
                                    tol = 1e-8,
                                    max_iter = 2000L,
                                    verbose = FALSE,
                                    max_k_fraction = 0.5,
                                    k_initial = 5L,
                                    eigs_opts = list()) {
  R <- symm(R)
  n <- nrow(R)
  T <- length(u_list)
  for (i in seq_len(T)) {
    if (length(u_list[[i]]) != n) stop("u vectors must all have length nrow(R)")
  }
  
  X <- R
  Y <- replicate(T, matrix(0, n, n), simplify = FALSE)
  iter <- 0L
  converged <- FALSE
  rel_changes <- numeric(0)
  dense_used_count <- 0L
  start_time <- proc.time()[3]
  
  repeat {
    if (iter >= max_iter) break
    X_old <- X
    for (t in seq_len(T)) {
      u <- as.numeric(u_list[[t]])
      V <- X - Y[[t]]
      pr <- project_psd_shifted_iterative(V, u,
                                          max_k_fraction = max_k_fraction,
                                          k_initial = k_initial,
                                          eigs_opts = eigs_opts)
      P <- pr$P
      if (isTRUE(pr$used_dense)) dense_used_count <- dense_used_count + 1L
      Y[[t]] <- P - V
      X <- P
    }
    iter <- iter + 1L
    diff <- X - X_old
    rel_change <- norm(diff, type = "F") / max(1.0, norm(X_old, type = "F"))
    rel_changes <- c(rel_changes, rel_change)
    if (verbose && (iter %% 5 == 0 || rel_change < tol)) {
      cat(sprintf("[iter %d] rel_change = %0.3e, dense_used_count = %d\n",
                  iter, rel_change, dense_used_count))
    }
    if (rel_change <= tol) {
      converged <- TRUE
      break
    }
  }
  
  time_sec <- proc.time()[3] - start_time
  info <- list(converged = converged,
               iterations = iter,
               time_sec = time_sec,
               rel_changes = rel_changes,
               dense_used_count = dense_used_count)
  return(list(S = symm(X), info = info))
}

# -------------------------
# Demo / quick test
# -------------------------

set.seed(2)
n <- 120
T <- 8

# Construct PSD R
A <- matrix(rnorm(n * floor(n/4)), n, floor(n/4))
R <- A %*% t(A)

# create rank-one W_t = u u^T (scaled random vectors)
u_list <- vector("list", length = T)
for (i in seq_len(T)) {
  v <- rnorm(n)
  v <- v / sqrt(sum(v^2))
  u_list[[i]] <- v * (0.5 + runif(1) * 1.5)
}

res <- dykstra_rank1_iterative(R, u_list, tol = 1e-6, max_iter = 100L, verbose = TRUE,
                               max_k_fraction = 0.6, k_initial = 4L)
cat("Converged?", res$info$converged, "iters:", res$info$iterations,
    "time (s):", res$info$time_sec, "\n")
cat("dense_used_count:", res$info$dense_used_count, "\n")

S_proj <- res$S
# Verify constraints S - W_t PSD (check minimum eigenvalue)
worst <- 0
viol <- FALSE
for (i in seq_len(T)) {
  M <- symm(S_proj - tcrossprod(u_list[[i]], u_list[[i]]))
  evs <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(evs)
  worst <- min(worst, min_eig)
  if (min_eig < -1e-7) viol <- TRUE
}
cat("Most negative eigenvalue across constraints:", worst, "\n")
cat("Any significant violation?", viol, "\n")
cat("||S - R||_F =", norm(S_proj - R, type = "F"), "\n")
