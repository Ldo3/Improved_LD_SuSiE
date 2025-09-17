library(susieR)
library(Matrix)
library(RSpectra)

symm <- function(A) {
  # enforce symmetry
  return((A + t(A)) / 2)
}

## one trait, full SVD
proj_B <- function(R, v){
  vvT = tcrossprod(v)
  R_minus = R - vvT
  eigen_R = eigen(R_minus, symmetric = TRUE)

  V <- eigen_R$vectors
  D_plus <- pmax(eigen_R$values, 0)

  R_plus <- V %*% (D_plus * t(V)) + vvT
  return(R_plus)
}

## one trait, full SVD, Dysktra (guarantee diagonal near 1)
proj_Dykstra <- function(R,
                         v,
                         max_iter=200L,
                         tol=1e-8){

  iter = 0L
  rel_changes <- numeric(0)
  start_time <- proc.time()[3]

  J = nrow(R)
  if (length(v) != J) stop("v vector must have length nrow(R)")
  Q <- replicate(2, matrix(0, J, J), simplify = FALSE) ## store minus residual in Dykstra algorithm

  repeat{
    if (iter >= max_iter) break

    R_old <- R
    ## Project diag
    R_1 <- R - Q[[1]]  ## x_k + p_k
    R_2 <- R - Q[[1]]
    diag(R_1) <- 1  ## y_k

    Q[[1]] <- R_1 - R_2  ## - p_{k+1}
    ## Project B
    R_1 = R_1 - Q[[2]]  ## y_k + q_k
    R = proj_B(R_1, v)  ## P(y_k + q_k)
    Q[[2]] <- R - R_1

    iter <- iter + 1L
    diff <- R - R_old
    rel_change <- norm(diff, type = "F") / max(1.0, norm(R_old, type = "F"))
    rel_changes <- c(rel_changes, rel_change)
    if (rel_change <= tol) {
      converged <- TRUE
      break
    }
  }

  time_sec <- proc.time()[3] - start_time
  info <- list(time_sec = time_sec,
               rel_changes = rel_changes)
  diag(R) <- 1
  return(list(R = symm(R), info = info))
}

## multi trait, full SVD, Dysktra (guarantee diagonal near 1)
proj_Dykstra_mtraits <- function(R,
                                 v_lists,
                                 max_iter=200L,
                                 tol=1e-8){
  iter = 0L
  rel_changes <- numeric(0)
  start_time <- proc.time()[3]

  J = nrow(R)
  Tr <- length(v_lists)
  for (t in 1:Tr){
    if (length(v_lists[[t]]) != J) stop("All v vectors must have length nrow(R)")
  }

  Q <- replicate(Tr, matrix(0, J, J), simplify = FALSE) ## store minus residual in Dykstra algorithm
  QD <- replicate(Tr, matrix(0, J, J), simplify = FALSE) ## for project Diag

  repeat{
    if (iter >= max_iter) break

    R_old <- R
    for (t in 1:Tr){
      ## Project diag
      R_1 <- R - QD[[t]]  ## x_k + p_k
      R_2 <- R - QD[[t]]
      diag(R_1) <- 1  ## y_k

      QD[[t]] <- R_1 - R_2  ## - p_{k+1}
      ## Project B
      v <- as.numeric(v_lists[[t]])
      R_1 <- R_1 - Q[[t]]  ## y_k + q_k
      R <- proj_B(R_1, v)  ## P(y_k + q_k)
      Q[[t]] <- R - R_1
    }

    iter <- iter + 1L
    diff <- R - R_old
    rel_change <- norm(diff, type = "F") / max(1.0, norm(R_old, type = "F"))
    rel_changes <- c(rel_changes, rel_change)
    if (rel_change <= tol) {
      converged <- TRUE
      break
    }
  }

  time_sec <- proc.time()[3] - start_time
  info <- list(time_sec = time_sec,
               rel_changes = rel_changes)
  return(list(R = symm(R), info = info))
}


###################################################################
# In the following, we use RSpectra instead of eigen
# to find only k smallest eigenvalue of R instead of all eigenvalues
# and truncate all negative eigenvalues among them
# If the k smallest eigenvalues do not contain all negative eigenvalues
# then increase k



## one trait, Dysktra, RSpectra
project_psd_shifted_iterative <- function(M, u,
                                          max_k_fraction = 0.5,
                                          k_initial = 5L,
                                          k_step_factor = 2L,
                                          tol_eig = 1e-4,
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
dykstra_rank1_iterative_opt1 <- function(R, v_list,
                                    tol = 1e-8,
                                    max_iter = 200L,
                                    verbose = FALSE,
                                    max_k_fraction = 0.5,
                                    k_initial = 5L,
                                    tol_eig = 1e-4,
                                    eigs_opts = list()) {
  R <- symm(R)
  J <- nrow(R)
  T <- length(v_list)
  for (i in seq_len(T)) {
    if (length(v_list[[i]]) != J) stop("v vectors must all have length nrow(R)")
  }

  Q <- replicate(T, matrix(0, J, J), simplify = FALSE)
  QD <- replicate(T, matrix(0, J, J), simplify = FALSE)
  iter <- 0L
  converged <- FALSE
  rel_changes <- numeric(0)
  dense_used_count <- 0L
  start_time <- proc.time()[3]

  repeat {
    if (iter >= max_iter) break
    R_old <- R

    for (t in seq_len(T)) {
      R1 <- R - QD[[t]]  ## x_k + p_k
      R2 <- R - QD[[t]]
      diag(R1) <- 1  ## y_k
      QD[[t]] <- R1 - R2  ## - p_{k+1}

      v <- as.numeric(v_list[[t]])
      R1 <- R1 - Q[[t]]
      pr <- project_psd_shifted_iterative(R1, v,
                                          max_k_fraction = max_k_fraction,
                                          tol_eig = tol_eig,
                                          k_initial = k_initial,
                                          eigs_opts = eigs_opts)
      R <- pr$P
      if (isTRUE(pr$used_dense)) dense_used_count <- dense_used_count + 1L
      Q[[t]] <- R - R1
    }
    iter <- iter + 1L
    diff <- R - R_old
    rel_change <- norm(diff, type = "F") / max(1.0, norm(R_old, type = "F"))
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
  return(list(R = symm(R), info = info))
}

########################################################

dykstra_rank1_iterative_opt2 <- function(R, v_list,
                                         tol = 1e-8,
                                         max_iter = 200L,
                                         k=10L) {
  R <- symm(R)
  J <- nrow(R)
  T <- length(v_list)
  for (i in seq_len(T)) {
    if (length(v_list[[i]]) != J) stop("v vectors must all have length nrow(R)")
  }

  Q <- replicate(T, matrix(0, J, J), simplify = FALSE)
  QD <- replicate(T, matrix(0, J, J), simplify = FALSE)
  iter <- 0L
  converged <- FALSE
  rel_changes <- numeric(0)
  start_time <- proc.time()[3]

  repeat {
    if (iter >= max_iter) break
    R_old <- R

    for (t in seq_len(T)) {
      R1 <- R - QD[[t]]  ## x_k + p_k
      R2 <- R - QD[[t]]
      diag(R1) <- 1  ## y_k
      QD[[t]] <- R1 - R2  ## - p_{k+1}

      R1 <- R1 - Q[[t]]
      v <- as.numeric(v_list[[t]])
      # Project R1
      A = R1 - tcrossprod(v)
      ev = eigs_sym(A, k = k, which = "SA")
      Pr = A - ev$vectors %*% (ev$values * t(ev$vectors)) ## remove k smallest eigenvalues
      R <- Pr + tcrossprod(v)
      # calculate residual
      Q[[t]] <- R - R1
    }
    iter <- iter + 1L
    diff <- R - R_old
    rel_change <- norm(diff, type = "F") / max(1.0, norm(R_old, type = "F"))
    rel_changes <- c(rel_changes, rel_change)
    if (rel_change <= tol) {
      converged <- TRUE
      break
    }
  }

  time_sec <- proc.time()[3] - start_time
  info <- list(converged = converged,
               iterations = iter,
               time_sec = time_sec,
               rel_changes = rel_changes)

  return(list(R = symm(R), info = info))
}






