library(susieR)
library(Matrix)

seed = 20
set.seed(seed)
setwd("~/Documents/Improved_LD_SuSiE")
# setwd("~/SuSiE-covariance/Improved_LD_SuSiE")

gtex = readRDS("data/Thyroid_ENSG00000132855.rds")
maf = apply(gtex, 2, function(x) sum(x)/2/length(x))
X0 = gtex[, maf > 0.01]
dim(X0)
X = na.omit(X0)
dim(X)
snp_total = ncol(X0)
n = nrow(X0)
p = 50
# Start from a random point on the genome
indx_start = sample(1: (snp_total - p), 1)
X = X0[, indx_start:(indx_start + p -1)]
# View(cor(X)[1:10, 1:10])

## sub-sample into two
out_sample = sample(1:n, 100)
X_out = X[out_sample, ]
X_in = X[setdiff(1:n, out_sample), ]
sum(is.na(X_out))

rm_p = c(which(diag(cov(X_in))==0), which(diag(cov(X_out))==0))
indx_p = setdiff(1:p, rm_p)
X_in = X_in[, indx_p]
X_out = X_out[, indx_p]
## out-sample LD matrix
Rp = cov(X_out)
R0 = cov(X_in)

library(ggplot2)
library(reshape2)

df1 <- melt(R0)
df2 <- melt(Rp)

p1 <- ggplot(df1, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  coord_fixed() +
  ggtitle("In-sample Covariance Matrix")

p2 <- ggplot(df2, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red") +
  coord_fixed() +
  ggtitle("Out-of-sample Covariance Matrix")

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)


#### log IW(R0 | nu0 * Rp, nu0 + J + 1)
log_multigamma_vec <- function(a, p) {
  # vectorized multivariate gamma
  j <- 1:p
  # sum over j, but broadcasting a over j
  (p*(p-1)/4)*log(pi) + 
    rowSums(matrix(lgamma(a), nrow=length(a), ncol=p, byrow=FALSE) + 
              matrix((1 - j)/2, nrow=length(a), ncol=p, byrow=TRUE))
}

log_iw <- function(R0, Rp, nu_vec) {
  p <- nrow(R0)
  R0 = R0 + 1e-8 * diag(rep(1, p))
  Rp = Rp + 1e-8 * diag(rep(1, p))
  # Precompute expensive shared quantities
  logdetRp <- determinant(Rp, logarithm = TRUE)$modulus + p * log(nu_vec)
  logdetR0   <- determinant(R0,   logarithm = TRUE)$modulus
  tr_term   <- nu_vec * sum(Rp * solve(R0))
  
  llhs = (.5 * (nu_vec + p + 1) * logdetRp
         - .5 * (nu_vec + p + 1) * p * log(2)
         - log_multigamma_vec((nu_vec + p + 1) / 2, p)
         - .5 * (nu_vec + p + 1) * logdetR0
         - .5 * tr_term)
  
  as.numeric(llhs)
}

nu_vec = c(1, 2, 3, 4)
llhs = log_iw(R0, Rp, nu_vec)









