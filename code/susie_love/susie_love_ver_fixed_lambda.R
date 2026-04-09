## wrapped R function for susie_love with fixed lambda
## Prior mean: R_bar = (1 - lambda) * R0 + lambda * diag(c)
## R0 can be substituted by R_tilde so that diag(R_tilde)
## R_tilde = diag(R)^{1/2} diag(R0)^{-1/2} R0 diag(R0)^{-1/2} diag(R)^{1/2}

# ---- 1. Overall loss function ----
## Loss function = - ELBO / J

## b_res = q(b) and prior sigma_0^2
## Z, nu_q, nu_0 = q(R) and prior parameter nu_0
## The mean R_q of q(R) is parametrized by Z to ensure diagonal constraint
## k = rank of Z

compute_loss <- function(b_res, Z, nu_q, nu_0, v, N, R_bar, c_diag = NULL) {
  # 1. Extract VI parameters
  alpha = b_res$alpha
  mu = b_res$mu
  sigma2_post = b_res$sigma2_post
  sigma02 = b_res$sigma02
  
  b_bar = alpha * mu
  b_bar_sum = colSums(b_bar)
  
  B = compute_E_bbT(alpha, mu, sigma2_post)
  
  J <- nrow(Z)
  k <- ncol(Z)
  
  if (is.null(c_diag)) c_diag <- rep(1, J)
  
  # ---------------------------------------------------------
  # 2. Reconstruct fast R_q components (FIXED VARIABLE NAME)
  # ---------------------------------------------------------
  s <- 1 + rowSums(Z^2)
  alpha_diag <- sqrt(c_diag / s) # Renamed to avoid collision
  
  F_mat <- Z * alpha_diag
  D_diag <- alpha_diag^2
  D_inv_diag <- 1 / D_diag
  
  F_T_Dinv_F <- t(F_mat) %*% (F_mat * D_inv_diag)
  M <- diag(k) + F_T_Dinv_F
  
  chol_M <- chol(M)
  log_det_M <- 2 * sum(log(diag(chol_M)))
  log_det_Rq <- sum(log(D_diag)) + log_det_M
  
  M_inv <- chol2inv(chol_M)
  D_inv_F <- F_mat * D_inv_diag
  Rq_inv <- diag(D_inv_diag) - D_inv_F %*% tcrossprod(M_inv, D_inv_F)
  
  # ---------------------------------------------------------
  # 3. Compute Pre-requisites and Traces
  # ---------------------------------------------------------
  log_det_Rbar <- as.numeric(determinant(R_bar, logarithm = TRUE)$modulus)
  V <- ((nu_q + J + 1) / nu_q) * (nu_0 * R_bar + N * tcrossprod(v))
  
  trace_Rq_B <- sum(crossprod(F_mat, B) * t(F_mat)) + sum(D_diag * diag(B))
  trace_Rq_inv_V <- sum(Rq_inv * V)
  
  # ---------------------------------------------------------
  # 4. Compute ELBO Terms
  # ---------------------------------------------------------
  term_1_and_2 <- diff_lmvgamma((nu_q + J + 1) / 2, (nu_0 + J + 1) / 2, J) / J
  term3 <- - ((nu_q - nu_0 - 1) / 2) * mvdigamma((nu_q + J + 1) / 2, J) / J
  
  term4 <- ((nu_0 + J + 1) / 2) * (log(nu_0) - log(nu_q))
  term5 <- ((nu_q + J + 1) / 2)
  term6 <- ((nu_0 + J + 1) / 2) * log_det_Rbar / J
  
  term_nu = term_1_and_2 + term3 + term4 + term5 + term6
  
  term7 <- - ((nu_0 + J + 2) / 2) * log_det_Rq 
  term8 <- - (N / 2) * trace_Rq_B 
  term9 <- - (1 / 2) * trace_Rq_inv_V 
  
  term_R = (term7 + term8 + term9) / J
  
  term_vb = N * v %*% b_bar_sum 
  
  # term KL(q(b) || p(b))
  alpha_safe <- pmax(alpha, 1e-12) # Safely references b_res$alpha
  kl_cat <- sum(alpha * (log(alpha_safe) - log(1 / J)))
  sigma02_mat <- matrix(sigma02, nrow = nrow(alpha), ncol = ncol(alpha))
  kl_norm <- 0.5 * sum(alpha * ( 
    (sigma2_post + mu^2) / sigma02_mat - 1 - log(sigma2_post / sigma02_mat) 
  ))
  term_kl_b <- kl_cat + kl_norm
  
  term_b = (term_vb - term_kl_b) / J
  
  elbo_val <- term_nu + term_R + term_b
  
  return(- elbo_val)
}


## In the following we extract the parts of loss function 
## with respect to R_q and nu to optimize them alternatively
## Update of q(b) exists in closed-form

# ---- 2. Subroutine to optimize R_q (parametrized by Z) ----

## the part of the loss function wrt R_q is
## a * log |R_q| + b * Tr(R_q B) + c * Tr(R_q^{-1} V)
## a = (nu_0 + J + 2) / (2 * J)
## b = N / (2 * J)
## c = 1 / (2 * J)


objective_Rq <- function(z_vec, B, V, c_diag, J, k, a_const, b_const, c_const) {
  
  # 1. Reconstruct Z, F, and Diagonal D
  Z <- matrix(z_vec, nrow = J, ncol = k)
  s <- 1 + rowSums(Z^2)
  alpha <- sqrt(c_diag / s)
  
  F_mat <- Z * alpha
  D_diag <- alpha^2
  D_inv_diag <- 1 / D_diag
  
  # 2. Woodbury / Sylvester speedups for R_q^{-1} and log|R_q|
  F_T_Dinv_F <- t(F_mat) %*% (F_mat * D_inv_diag)
  M <- diag(k) + F_T_Dinv_F
  
  chol_M <- chol(M)
  log_det_M <- 2 * sum(log(diag(chol_M)))
  log_det_Rq <- sum(log(D_diag)) + log_det_M
  
  M_inv <- chol2inv(chol_M)
  D_inv_F <- F_mat * D_inv_diag
  
  # Compute dense R_q^{-1} (Needed for trace(R_q^{-1} V) anyway)
  Rq_inv <- diag(D_inv_diag) - D_inv_F %*% tcrossprod(M_inv, D_inv_F)
  
  # 3. Compute Traces
  # trace(R_q B) = trace(F F^T B) + trace(D B)
  trace_Rq_B <- sum(crossprod(F_mat, B) * t(F_mat)) + sum(D_diag * diag(B))
  
  # trace(R_q^{-1} V)
  trace_Rq_inv_V <- sum(Rq_inv * V)
  
  # 4. Final Loss Value
  loss <- (a_const * log_det_Rq) + (b_const * trace_Rq_B) + (c_const * trace_Rq_inv_V)
  return(loss)
}

# --- 2. THE ANALYTICAL GRADIENT ---
gradient_Rq <- function(z_vec, B, V, c_diag, J, k, a_const, b_const, c_const) {
  
  # 1. Reconstruct parameters
  Z <- matrix(z_vec, nrow = J, ncol = k)
  s <- 1 + rowSums(Z^2)
  alpha <- sqrt(c_diag / s)
  
  F_mat <- Z * alpha
  D_diag <- alpha^2
  D_inv_diag <- 1 / D_diag
  
  # 2. Reconstruct R_q^{-1}
  F_T_Dinv_F <- t(F_mat) %*% (F_mat * D_inv_diag)
  M <- diag(k) + F_T_Dinv_F
  M_inv <- chol2inv(chol(M))
  
  D_inv_F <- F_mat * D_inv_diag
  Rq_inv <- diag(D_inv_diag) - D_inv_F %*% tcrossprod(M_inv, D_inv_F)
  
  # 3. Construct Matrix Gradient G
  # G = a * R_q^{-1} + b * B - c * R_q^{-1} V R_q^{-1}
  Rq_inv_V <- Rq_inv %*% V
  Rq_inv_V_Rq_inv <- Rq_inv_V %*% Rq_inv
  
  G <- (a_const * Rq_inv) + (b_const * B) - (c_const * Rq_inv_V_Rq_inv)
  G <- (G + t(G)) / 2 # Ensure strict numerical symmetry
  
  # 4. Compute gradient w.r.t Z
  U <- 2 * G %*% F_mat
  G_diag <- diag(G)
  
  dot_vals <- rowSums(U * Z) 
  scalar_term <- (alpha * dot_vals + 2 * G_diag * D_diag) / s
  
  grad_Z <- (U * alpha) - (Z * scalar_term)
  
  return(as.vector(grad_Z))
}

initialize_Z_from_Rq <- function(Rq_init, c_diag, k) {
  J <- nrow(Rq_init)
  
  # 1. Get rank-k approximation via Eigendecomposition
  eigen_res <- eigen(Rq_init, symmetric = TRUE)
  
  # Extract top k eigenvalues and eigenvectors
  evals <- eigen_res$values[1:k]
  evecs <- eigen_res$vectors[, 1:k, drop = FALSE]
  
  # Ensure eigenvalues are positive (floor them at a small number if not)
  evals <- pmax(evals, 1e-8)
  
  # Construct initial F matrix: V * sqrt(Lambda)
  F_init <- evecs %*% diag(sqrt(evals), nrow = k, ncol = k)
  
  # 2. Compute row norms squared of F
  F_row_norms_sq <- rowSums(F_init^2)
  
  # 3. Safety Catch: Ensure ||F_{i*}||^2 is strictly less than c_i
  # We cap it at 95% of c_i to ensure alpha remains strictly positive
  # and Z doesn't initialize with massive numbers.
  max_allowed <- 0.95 * c_diag
  
  for (i in 1:J) {
    if (F_row_norms_sq[i] > max_allowed[i]) {
      # Scale the row down so its squared norm equals max_allowed
      scale_factor <- sqrt(max_allowed[i] / F_row_norms_sq[i])
      F_init[i, ] <- F_init[i, ] * scale_factor
      F_row_norms_sq[i] <- max_allowed[i] # Update the norm
    }
  }
  
  # 4. Reverse the mapping to solve for alpha and Z
  alpha <- sqrt(c_diag - F_row_norms_sq)
  Z_init <- F_init / alpha  # Divides row-wise automatically in R
  
  # Return as the flattened vector expected by optim()
  return(as.vector(Z_init))
}

# ---- 3. Objective function of nu_q and nu_0 ----
objective_nu <- function(nu_vec, J, N, log_det_Rbar, log_det_Rq, trace_Rq_inv_Rbar, trace_Rq_inv_vvT) {
  nu_q <- nu_vec[1]
  nu_0 <- nu_vec[2]
  
  # 1. Terms from term_nu
  term_1_and_2 <- diff_lmvgamma((nu_q + J + 1) / 2, (nu_0 + J + 1) / 2, J) / J
  term3 <- - ((nu_q - nu_0 - 1) / 2) * mvdigamma((nu_q + J + 1) / 2, J) / J
  term4 <- ((nu_0 + J + 1) / 2) * (log(nu_0) - log(nu_q))
  term5 <- ((nu_q + J + 1) / 2)
  term6 <- ((nu_0 + J + 1) / 2) * log_det_Rbar / J
  
  term_nu <- term_1_and_2 + term3 + term4 + term5 + term6
  
  # 2. Terms from term_R
  term7 <- - ((nu_0 + J + 2) / 2) * log_det_Rq
  
  # trace(Rq_inv %*% V) expanded
  trace_V_base <- nu_0 * trace_Rq_inv_Rbar + N * trace_Rq_inv_vvT
  term9 <- - (1 / 2) * ((nu_q + J + 1) / nu_q) * trace_V_base
  
  term_R <- (term7 + term9) / J
  
  # Return the negative ELBO contribution
  return(-(term_nu + term_R))
}


# ---- 4. Supporting functions ----

diff_lmvgamma <- function(a_q, a_0, J) {
  # Computes log Gamma_J(a_q) - log Gamma_J(a_0)
  # The (J(J-1)/4)*log(pi) constant is analytically canceled.
  
  # Create the shift sequence once: (1-i)/2
  shift_seq <- (1 - 1:J) / 2
  
  # Compute the sum of the differences of the 1D log-gamma functions
  # Vectorized subtraction is much faster and more numerically stable
  return(sum(lgamma(a_q + shift_seq) - lgamma(a_0 + shift_seq)))
}

mvdigamma <- function(a, p) {
  # psi_p(a) = sum_{i=1}^p psi(a + (1-i)/2)
  # Note: R's built-in digamma() represents the psi function
  sum_terms <- sum(digamma(a + (1 - 1:p) / 2))
  return(sum_terms)
}

compute_E_bbT <- function(alpha, mu, sigma2_post) {
  L <- nrow(alpha)
  J <- ncol(alpha)
  
  # Expected effect vectors for each L (L x J matrix)
  b_bar_l <- alpha * mu
  
  # Total expected effect vector (J x 1 vector)
  b_bar <- colSums(b_bar_l)
  
  # 1. Total outer product: b_bar %*% b_bar^T
  b_bar_b_bar_T <- tcrossprod(b_bar) 
  
  # Handle sigma2 formatting (it can be a J-length vector if fixed across L, 
  # or an L x J matrix if prior variance V differs per L)
  
  # 2. Diagonal matrix D: sum of second moments
  D_vec <- colSums(alpha * (mu^2 + sigma2_post))
  D <- diag(D_vec)
  
  # 3. Sum of the within-component outer products
  sum_bl_bl_T <- matrix(0, nrow = J, ncol = J)
  for (l in 1:L) {
    sum_bl_bl_T <- sum_bl_bl_T + tcrossprod(b_bar_l[l, ])
  }
  
  # Combine to get the final expectation E_q[bb^T]
  E_bbT <- D + b_bar_b_bar_T - sum_bl_bl_T
  
  return(E_bbT)
}

# ---- 5. Initialization Z ----

## Initialize Z by minimizing KL(IW(R_q), IW(R_bar))
## R_bar is the prior mean 

initialize_Z_from_prior_R <- function(R_bar, nu_0, nu_q, c_diag, k) {
  
  ## 1. Get initial estimate of Z via SVD
  J <- nrow(R_bar)
  
  # (a). Get rank-k approximation via Eigendecomposition
  eigen_res <- eigen(R_bar, symmetric = TRUE)
  
  # Extract top k eigenvalues and eigenvectors
  evals <- eigen_res$values[1:k]
  evecs <- eigen_res$vectors[, 1:k, drop = FALSE]
  
  # Ensure eigenvalues are positive (floor them at a small number if not)
  evals <- pmax(evals, 1e-8)
  
  # Construct initial F matrix: V * sqrt(Lambda)
  F_init <- evecs %*% diag(sqrt(evals), nrow = k, ncol = k)
  
  # (b). Compute row norms squared of F
  F_row_norms_sq <- rowSums(F_init^2)
  
  # (c). Safety Catch: Ensure ||F_{i*}||^2 is strictly less than c_i
  # We cap it at 95% of c_i to ensure alpha remains strictly positive
  # and Z doesn't initialize with massive numbers.
  max_allowed <- 0.95 * c_diag
  
  for (i in 1:J) {
    if (F_row_norms_sq[i] > max_allowed[i]) {
      # Scale the row down so its squared norm equals max_allowed
      scale_factor <- sqrt(max_allowed[i] / F_row_norms_sq[i])
      F_init[i, ] <- F_init[i, ] * scale_factor
      F_row_norms_sq[i] <- max_allowed[i] # Update the norm
    }
  }
  
  # (d). Reverse the mapping to solve for alpha and Z
  alpha <- sqrt(c_diag - F_row_norms_sq)
  Z_init <- F_init / alpha  # Divides row-wise automatically in R
  
  # Return as the flattened vector expected by optim()
  z_init = as.vector(Z_init)
  # 2. Get a better initialization by minimizing the KL 
  ## Note that we can code a separate function for KL
  ## but I am so lazy so I just use the objective function 
  ## of Rq and set v = 0, B = 0 (which becomes the KL)
  
  scalar_V <- (nu_q + J + 1) / nu_q
  V_matrix <- scalar_V * (nu_0 * R_bar)
  a_const <- (nu_0 + J + 2) / (2 * J)
  b_const <- N / (2 * J)
  c_const <- 1 / (2 * J)
  B_0 = matrix(0, nrow = J, ncol = J)  ## 0
  v_0 = rep(0, J) ## 0
  
  optim_res <- optim(
    par = z_init,
    fn = objective_Rq,
    gr = gradient_Rq,
    B = B_0,          
    V = V_matrix,           
    c_diag = c_diag, J = J, k = k,
    a_const = a_const, b_const = b_const, c_const = c_const,
    method = "L-BFGS-B",
    control = list(maxit = 30) 
  )
  
  z_init <- optim_res$par
  
  return(z_init)
}


