SER <- function(V_xx, V_xy, n, sigma2, sigma02){
  # Single Effect Regression (SER) model
  # Note that we already centered and standardize X and y so that x_j^T x_j = y^T y = 1
  # gamma ~ Multi(1, pi); pi uniform
  # y | gamma=1 ~ N(X_j*b_j, sigma2) 
  # b_j ~ N(0, sigma02) 
  # Summarize in V_xx (p x p) cov matrix and V_xy (p) inner product between X_j and y.
  # 
  # Return: alpha (posterior inclusion probability; PIP), 
  # mus (posterior mean of b), 
  # sigmas (posterior variance of b)
  
  b_hat = V_xy
  s2 = sigma2 / n
  sigma12 = s2 * sigma02 / (s2 + sigma02)
  mus = sigma12 / s2 * b_hat 
  
  log_BFs = ((log(s2) - log(s2 + sigma02)) / 2 
             + b_hat^2 / (2*s2) * sigma02 / (s2 + sigma02))
  
  unnorm_logBFs = log_BFs - max(log_BFs)
  unnorm_BFs = exp(unnorm_logBFs)
  alpha = unnorm_BFs / sum(unnorm_BFs)
  
  return(list(alpha=alpha, mus=mus, sigma12=sigma12))
}

SuSiE_rss <- function(V_xx, X_xy, n, L, max_iter, sigma02, sigma2){
  b_bar = matrix(0, nrow = L, ncol = J)
  b_bar2 = matrix(0, nrow = L, ncol = J)
  alphas = matrix(0, nrow = L, ncol = J)
  mus = matrix(0, nrow = L, ncol = J)
  sigmas = matrix(0, nrow = L, ncol = J)
  
  for (iter in (1:max_iter)){
    V_xy_bar = V_xy - X_xx %*% colSums(b_bar)  ## residual signal
    
    for (ell in 1:L){
      V_xy_bar_ell = V_xy_bar + V_xx %*% b_bar[ell, ]  ## add back ell-th signal
      ret = SER(V_xx, V_xy_bar_ell, n, sigma2, sigma02)
      alphas[ell, ] = ret$alpha
      mus[ell, ] = ret$mus
      sigmas[ell, ] = ret$sigma12
      b_bar[ell, ] = alpha_ell * mus_ell
      b_bar2[ell, ] = alpha_ell * (mus_ell^2 + sigmas_ell)
      V_xy_bar = V_xy_bar_ell - V_xx %*% b_bar[ell, ]
    }
    
    
  }
  return(list(alphas = alphas, mus = mus, sigmas = sigmas))
}


