A <- tcrossprod(zhat)
A <- A/max(abs(A))
w <- 1
# w <- 0.5
R2 <- w*R_out + (1-w)*A
fit_out <- susie_rss(zhat,R2,L = 4,n = n,estimate_prior_method = "EM",
                     tol = 1e-6,min_abs_corr = 0,verbose = TRUE)
