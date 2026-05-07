library(MCMCpack)

# Simulate a data set from this model:
#
#   x ~ N(0,S)
#   S ~ IW(v,S0)
set.seed(1)
S0 <- rbind(c(1,0.95,-0.1),
            c(0.95,1,0),
            c(-0.1,0,1))
