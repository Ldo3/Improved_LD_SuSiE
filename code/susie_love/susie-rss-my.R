compute_elbo <- function(v, R, N, alpha, mu, sigma2_post, sigma02) {
  L <- nrow(alpha)
  J <- ncol(alpha)
  
  # 1. Expected Log-Likelihood
  b_bar_l <- alpha * mu          # Expected effects per component (L x J)
  b_bar <- colSums(b_bar_l)      # Total expected effect (J)
  
  # Calculate tr(R * E[bb^T])
  D_vec <- colSums(alpha * (mu^2 + sigma2_post)) 
  
  term_D <- sum(diag(R) * D_vec)
  term_cross_total <- sum(b_bar * (R %*% b_bar))
  term_cross_within <- sum(b_bar_l * t(R %*% t(b_bar_l))) # Vectorized sum(b_l^T R b_l)
  
  expected_ll <- sum(N * v * b_bar) - 0.5 * N * (term_D + term_cross_total - term_cross_within)
  
  # 2. KL Divergence sum_l KL( q(b_l) || p(b_l) )
  # Vectorize the Categorical KL divergence
  alpha_safe <- pmax(alpha, 1e-12)
  kl_cat <- sum(alpha * log(alpha_safe * J))
  
  # Replicate the length-L sigma02 vector into an L x J matrix for element-wise math
  sigma02_mat <- matrix(sigma02, nrow = nrow(alpha), ncol = ncol(alpha))
  
  # Vectorize the Normal KL divergence
  kl_norm <- 0.5 * sum(alpha * ( 
    (sigma2_post + mu^2) / sigma02_mat - 1 - log(sigma2_post / sigma02_mat) 
  ))
  
  # Total KL divergence
  kl_div <- kl_cat + kl_norm
  
  # ELBO = Expected Log-Likelihood - KL Divergence
  return(expected_ll - kl_div)
}

susie_rss_my <- function(b_res, v, R, N, max_iter = 100, tol = 1e-4) {
  J <- length(v)
  R_diag <- diag(R)
  
  # Initialize variational parameters from the passed list
  alpha = b_res$alpha
  mu = b_res$mu
  sigma2_post = b_res$sigma2_post
  sigma02 = b_res$sigma02
  
  # Explicitly define L so the function doesn't rely on the global environment
  L <- nrow(alpha) 
  
  b_bar <- alpha * mu 
  elbo_history <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    for (l in 1:L) {
      if (L == 1) {
        b_other <- rep(0, J)
      } else {
        b_other <- colSums(b_bar[-l, , drop = FALSE])
      }
      
      r_l <- as.vector(v - (R %*% b_other))
      sigma2_post_l = 1 / (N * R_diag + 1 / sigma02[l])
      mu_l <- sigma2_post_l * N * r_l
      log_bf <- 0.5 * log(sigma2_post_l / sigma02[l]) + 0.5 * (mu_l^2 / sigma2_post_l)
      
      max_log_bf <- max(log_bf)
      w <- exp(log_bf - max_log_bf)
      alpha[l, ] <- w / sum(w)
      
      mu[l, ] <- mu_l
      b_bar[l, ] <- alpha[l, ] * mu[l, ]
      sigma2_post[l, ] = sigma2_post_l
      sigma02[l] = max(1e-8, sum(alpha[l, ] * (mu[l, ]^2 + sigma2_post[l, ])))
    }
    
    current_elbo <- compute_elbo(v, R, N, alpha, mu, sigma2_post, sigma02)
    elbo_history[iter] <- current_elbo
    
    ## check convergence using ELBO, uncomment if want to use
    # if (iter > 1) {
    #   elbo_diff <- current_elbo - elbo_history[iter - 1]
    #   if (elbo_diff < -1e-6) {
    #     warning(sprintf("Iteration %d: ELBO decreased by %g", iter, elbo_diff))
    #   }
    #   if (abs(elbo_diff) < tol) {
    #     cat(sprintf("Converged at iteration %d. Final ELBO: %.2f\n", iter, current_elbo))
    #     elbo_history <- elbo_history[1:iter]
    #     break
    #   }
    # }
  }
  
  pip <- 1 - apply(1 - alpha, 2, prod)
  
  return(list(
    alpha = alpha,
    mu = mu,
    sigma2_post = sigma2_post,
    sigma02 = sigma02,
    pip = pip,
    elbo = elbo_history
  ))
}


get_cs <- function(alpha, R, coverage = 0.95, min_abs_corr = 0.5) {
  L <- nrow(alpha)
  J <- ncol(alpha)
  
  # FIX: Convert the covariance matrix to a correlation matrix
  R_corr <- cov2cor(R) 
  
  cs_list <- list()
  cs_purity <- numeric()
  
  for (l in 1:L) {
    # Sort the PIPs for component l in descending order
    pip_l <- alpha[l, ]
    sorted_indices <- order(pip_l, decreasing = TRUE)
    sorted_pips <- pip_l[sorted_indices]
    
    # Find the minimum number of variables needed to reach the coverage threshold
    cum_pips <- cumsum(sorted_pips)
    cs_size <- which(cum_pips >= coverage)[1]
    
    
    if (is.na(cs_size)) next 
    cs_vars <- sorted_indices[1:cs_size]
    
    if (length(cs_vars) == 1) {
      purity <- 1.0 
    } else {
      # FIX: Use the correlation matrix subset
      R_sub <- R_corr[cs_vars, cs_vars] 
      purity <- min(abs(R_sub))
    }
    
    if (purity >= min_abs_corr) {
      cs_name <- paste0("L", l)
      cs_list[[cs_name]] <- cs_vars
      cs_purity[cs_name] <- purity
    }
  }
  return(list(cs = cs_list, purity = cs_purity))
}

plot_cs <- function(pip, cs_obj, true_indices = NULL, 
                    title = "SuSiE Credible Sets with Causal Ground Truth") {
  J <- length(pip)
  
  # 1. Base Plot: All SNPs as neutral grey dots (slightly smaller)
  # Keeping ylim slightly expanded to 1.15 to prevent clipping at the very top
  plot(1:J, pip, 
       pch = 16,          # Solid circles 
       col = "grey70",    # Neutral background points
       cex = 0.8,         
       ylim = c(0, 1.15), 
       ylab = "Posterior Inclusion Probability (PIP)", 
       xlab = "Variable Index",
       main = title,
       bty = "l",         # L-shaped bounding box (standard SuSiE look)
       las = 1)           # Horizontal y-axis labels
  
  # Add a grid for readability
  grid(col = "grey90")
  
  # 2. Add open circles around the TRUE causal SNPs (if provided)
  if (!is.null(true_indices)) {
    # Using open circles (pch = 1) that are larger (cex = 1.8) 
    # and bolder (lwd = 2) to encase the filled dots
    points(true_indices, pip[true_indices], 
           pch = 1,          # Open circle
           col = "black", 
           cex = 1.8,        # Large size to encase other points
           lwd = 2)          # Thick border
  }
  
  # 3. Layer the learned Credible Sets using the classic SuSiE palette
  cs_colors <- c("dodgerblue", "firebrick", "forestgreen", "darkorange", "purple", "cyan", "magenta")
  cs_list <- cs_obj$cs
  
  if (length(cs_list) > 0) {
    for (i in seq_along(cs_list)) {
      cs_vars <- cs_list[[i]]
      col_idx <- ((i - 1) %% length(cs_colors)) + 1 
      
      # Overplot CS variables as larger, brightly colored, solid points (pch = 16)
      points(cs_vars, pip[cs_vars], 
             col = cs_colors[col_idx], 
             pch = 16,       # Solid circles
             cex = 1.2)      
    }
    
    # Unobtrusive legend for Credible Sets
    legend("topright", 
           legend = paste("CS", names(cs_list)), 
           col = cs_colors[1:length(cs_list)], 
           pch = 16, 
           pt.cex = 1.2,
           bty = "n",
           cex = 0.85)    # Slightly smaller text
  } else {
    legend("topright", legend = "No pure CS found", bty = "n", cex = 0.85)
  }
  
  # 4. Legend for Causal Ground Truth
  if (!is.null(true_indices)) {
    legend("topleft", 
           legend = "True Causal SNP", 
           col = "black", 
           pch = 1,          # Open circle matching the plot
           pt.lwd = 2,       # Thick border matching the plot
           pt.cex = 1.5,     
           bty = "n",
           cex = 0.85)
  }
}



susie_plot_alpha <- function(alpha, R = NULL, b = NULL, main = "SuSiE PIP Visualization", 
                             coverage = 0.95, min_abs_corr = 0.5, ...) {
  
  if (!is.matrix(alpha)) {
    stop("alpha must be an L x p matrix")
  }
  
  # 1. Calculate overall PIPs
  p <- 1 - apply(1 - alpha, 2, prod)
  
  # 2. Identify Credible Sets (CS) and calculate purity
  cs_list <- list()
  
  for (l in 1:nrow(alpha)) {
    alpha_l <- alpha[l, ]
    ord <- order(alpha_l, decreasing = TRUE)
    cum_alpha <- cumsum(alpha_l[ord])
    idx <- which(cum_alpha >= coverage)[1]
    
    if (max(alpha_l) > 1e-4) {
      cs_vars <- ord[1:idx]
      
      # If R is provided, we still check purity for filtering (SuSiE standard)
      if (!is.null(R)) {
        if (length(cs_vars) == 1) {
          purity <- 1.0
        } else {
          R_sub <- abs(R[cs_vars, cs_vars])
          purity <- min(R_sub[upper.tri(R_sub)])
        }
        if (purity >= min_abs_corr) {
          cs_list[[paste0("L", l)]] <- cs_vars
        }
      } else {
        cs_list[[paste0("L", l)]] <- cs_vars
      }
    }
  }
  
  # Setup base plot parameters
  pos <- 1:length(p)
  if (is.null(b)) {
    b <- rep(0, length(p))
  }
  
  args <- list(...)
  if (!exists("xlab", args)) args$xlab <- "variable"
  if (!exists("ylab", args)) args$ylab <- "PIP"
  if (!exists("pch", args)) args$pch <- 16
  args$x <- pos
  args$y <- p
  args$main <- main
  
  # 3. Draw base PIP plot
  do.call(plot, args)
  
  # Original colors
  color <- c(
    "dodgerblue2", "green4", "#6A3D9A", "#FF7F00", "gold1",
    "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F",
    "gray70", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", 
    "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  
  legend_text <- list(col = vector(), cs_label = vector())
  
  # 4. Overlay Credible Set Colors
  for (i in rev(seq_along(cs_list))) {
    cs_vars <- cs_list[[i]]
    cs_name <- names(cs_list)[i]
    
    if (length(cs_vars) > 0) {
      col_to_use <- head(color, 1)
      points(pos[cs_vars], p[cs_vars], col = col_to_use, cex = 1.5, lwd = 2.5)
      
      legend_text$col <- append(col_to_use, legend_text$col)
      legend_text$cs_label <- append(cs_name, legend_text$cs_label)
      
      color <- c(color[-1], color[1])
    }
  }
  
  # 5. Simplified Legend (No size or purity displayed)
  if (length(legend_text$col) > 0) {
    legend("topright", 
           legend = legend_text$cs_label, 
           bty = "n", 
           col = legend_text$col, 
           cex = 0.65, 
           pch = 15)
  }
  
  # 6. Highlight True Effects
  points(pos[b != 0], p[b != 0], col = 2, pch = 16)
  
  return(invisible(cs_list))
}