source("SuSiE_rss.R")
n = 2000
J = 200

set.seed(1)
X = matrix(rnorm(n*J), nrow = n, ncol = J)
gamma = sample(c(1:J), size = 1)

b = rnorm(1) * 3

y = X[ , gamma] * b + rnorm(n)
print(gamma)

## standardize and compute sufficient statistics
X = scale(X) / sqrt(n-1)
y = scale(y) / sqrt(n-1)

V_xx = t(X) %*% X
V_xy = t(X) %*% y

ret = SER(V_xx = V_xx, V_xy = V_xy, n = n, sigma2 = 1, sigma02 = sqrt(3))
alpha = ret$alpha
print(alpha[gamma])
print(which.max(alpha))
print(gamma)
print(b)
