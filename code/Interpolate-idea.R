library(susieR)
library(Matrix)    


# setwd("~/Documents/Improved_LD_SuSiE")
setwd("~/SuSiE-covariance/Improved_LD_SuSiE")

source("code/R_algorithms.R")
gtex = readRDS("data/Thyroid_ENSG00000132855.rds")

num_causal_SNP = 1
num_reps = 100

num_covs = 12
coverages = matrix(0, nrow = num_reps, ncol = num_covs)
powers = matrix(0, nrow = num_reps, ncol = num_covs)
no_CSs = matrix(0, nrow = num_reps, ncol = num_covs)
no_SNPs_CS = matrix(0, nrow = num_reps, ncol = num_covs)
compare_R = rep(0, num_reps)


# Remove SNPs with MAF < 0.01
maf = apply(gtex, 2, function(x) sum(x)/2/length(x))
X0 = gtex[, maf > 0.01]
dim(X0)
X = na.omit(X0)
dim(X)
snp_total = ncol(X0)

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
  true_effect = rnorm(num_causal_SNP) * 3
  beta[true_causal_SNPs] <- true_effect
  # plot(beta, pch=16, ylab='effect size')
  y <- X_in %*% beta + rnorm(n)
  y = scale(y)
  
  ## 0. compute summary statistics 
  sumstats <- univariate_regression(X_in, y)
  z_scores <- sumstats$betahat / sumstats$sebetahat
  # susie_plot(z_scores, y = "z", b=beta)
  min_cor = 0.5
  ## 1. fit the susie-rss model with in-sample R
  fitted_rss1 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = R, var_y = var(y), L = 10,
                           estimate_residual_variance = TRUE,
                           min_abs_corr=min_cor)
  # summary(fitted_rss1)$cs
  # p1 = susie_plot(fitted_rss1, y="PIP", b=beta)
  
  ## 2. fit the model with out-sample R
  fitted_rss2 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = R_hat, var_y = var(y), L = 10,
                           estimate_residual_variance = FALSE,
                           min_abs_corr=min_cor)
  
  
  ## 3. fit the model with out-sample R and L = 1
  fitted_rss3 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = R_hat, var_y = var(y), L = 1,
                           estimate_residual_variance = FALSE,
                           min_abs_corr=min_cor)
  lambda = 0.1
  fitted_rss4 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = ((1-lambda)*R_hat + lambda*diag(1, J)), var_y = var(y), L = 10,
                           estimate_residual_variance = FALSE,
                           min_abs_corr=min_cor)
  
  ## 3. P_1(R_hat)
  ## 6. P_1(R_hat)
  v = (t(X_in) %*% y) / (n-1)
  
  ret = proj_Dykstra(R=R_hat, v=v)
  R_hat_plus_diag1 = ret$R
  info1 = ret$info
  
  zero_mat = matrix(0, nrow=J, ncol=J)
  ret = proj_Dykstra(R=zero_mat, v=v)
  proj_zero_1 = ret$R
  
  fitted_rss5 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = R_hat_plus_diag1, var_y = var(y), L = 10,
                           estimate_residual_variance = F,
                           min_abs_corr=min_cor)
  
  fitted_rss6 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = proj_zero_1, var_y = var(y), L = 10,
                           estimate_residual_variance = F,
                           min_abs_corr=min_cor)
  
  fitted_rss7 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = 0.5*R_hat_plus_diag1 + 0.5*proj_zero_1, var_y = var(y), L = 10,
                           estimate_residual_variance = F,
                           min_abs_corr=min_cor)
  
  ## 4. P_5(R_hat)
  ## 7. P_5(0)
  
  v_lists = list(v)
  for (t in 1:9){
    beta_ = rep(0, J)
    causal_SNPs = sample(1:J, size=num_causal_SNP, replace=FALSE)
    true_effect = rnorm(num_causal_SNP) * 3
    beta_[causal_SNPs] <- true_effect
    # plot(beta, pch=16, ylab='effect size')
    y_ <- X_in %*% beta_ + rnorm(n)
    y_ = scale(y_)
    v_ = (t(X_in) %*% y_) / (n-1)
    v_lists <- c(v_lists, list(v_))
  }
  
  ret = proj_Dykstra_mtraits(R_hat, v_lists, max_iter=200L, tol=1e-6)
  R_hat_10 = ret$R
  info10 = ret$info
  
  ret = proj_Dykstra_mtraits(zero_mat, v_lists, max_iter=200L, tol=1e-6)
  proj_zero_10 = ret$R
  
  fitted_rss8 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = R_hat_10, var_y = var(y), L = 10,
                           estimate_residual_variance = F,
                           min_abs_corr=min_cor)
  fitted_rss9 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                           R = proj_zero_10, var_y = var(y), L = 10,
                           estimate_residual_variance = F,
                           min_abs_corr=min_cor)
  
  fitted_rss10 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                            R = 0.5*R_hat_10 + 0.5*proj_zero_10, var_y = var(y), L = 10,
                            estimate_residual_variance = F,
                            min_abs_corr=min_cor)
  
  fitted_rss11 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                            R = 0.75*R_hat_10 + 0.25*proj_zero_1, var_y = var(y), L = 10,
                            estimate_residual_variance = F,
                            min_abs_corr=min_cor)
  
  fitted_rss12 <- susie_rss(bhat = sumstats$betahat, shat = sumstats$sebetahat, n = n, 
                            R = 0.25*R_hat_10 + 0.75*proj_zero_1, var_y = var(y), L = 10,
                            estimate_residual_variance = F,
                            min_abs_corr=min_cor)
  
  # summary(fitted_rss5)$cs
  # susie_plot(fitted_rss5, y="PIP", b=beta)
  
  
  
  fitted_rss = list(fitted_rss1, 
                    fitted_rss2, 
                    fitted_rss3, 
                    fitted_rss4, 
                    fitted_rss5, 
                    fitted_rss6,
                    fitted_rss7, 
                    fitted_rss8,
                    fitted_rss9,
                    fitted_rss10,
                    fitted_rss11,
                    fitted_rss12
  )
  
  for (v in 1:num_covs){
    L_infer = nrow(summary(fitted_rss[[v]])$cs)
    if (is.null(summary(fitted_rss[[v]])$cs)) {
      coverages[seed, v] = -1
      powers[seed, v] = -1
      no_CSs[seed, v] = -1
      no_SNPs_CS[seed, v] = -1
    } else{
      no_contains = 0
      all_selected_SNPs = c()
      for (ell in 1:L_infer){
        this_CS = unlist(strsplit((summary(fitted_rss[[v]])$cs$variable[ell]), ",\\s*"))
        all_selected_SNPs = c(all_selected_SNPs, this_CS)
        no_contains = no_contains + (length(intersect(this_CS, true_causal_SNPs)) > 0)
        ## coverage = proportion of CS that contains a true causal SNP
        coverages[seed, v] = no_contains / L_infer
        ## power = proportion of causal SNP that is contained in a CS
        selected = length(intersect(all_selected_SNPs, true_causal_SNPs)) 
        powers[seed, v] = selected / num_causal_SNP
        
        ## number of CSs
        no_CSs[seed, v] = L_infer
        
        ## number of SNPs per CSs
        no_SNPs_CS[seed, v] = length(all_selected_SNPs) / L_infer
      }
    }
  }
}

colMeans(coverages)

library(ggplot2)

list_mat = list(coverages, powers, no_CSs, no_SNPs_CS)
list_name = c("Coverages", "Power", "Number of CSs", "Number of SNPs per CS")
plots = list()
for (i in 1:4){
  m = list_mat[[i]]
  m[m==-1] <- 0
  colnames(m) <- c('In-', 'Out-', 'L1', 'reg.', 
                   'P1(R)', 'P1(0)', 'P1(R/2)',
                   'P10(R)', 'P10(0)', 'P10(R/2)',
                   'int25', 'int75')
  means <- colMeans(m)
  sds   <- apply(m, 2, sd)
  df <- data.frame(
    variable = factor(colnames(m), levels = colnames(m)),
    mean = means,
    sd = sds
  )
  
  plots[[i]] = ggplot(df, aes(x = variable, y = mean)) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    labs(title = list_name[i],
         x = "Methods", y = "Mean ± SD")
}

library(patchwork)
wrap_plots(plots, ncol = 2)


View(no_CSs)


save(coverages, file='output/coverages_mat_multi_trait_L1.RData')
save(powers, file='output/powers_mat_multi_trait_L1.RData')
save(no_CSs, file='output/number_CSs_mat_multi_trait_L1.RData')
save(no_SNPs_CS, file='output/number_SNPs_per_CS_mat_multi_trait_L1.RData')





