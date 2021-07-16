setwd("/Users/xinyi/Documents/GitHub/network_ngc/R")
library('microbenchmark')
source("defn.R")
source("simulate.r")
source("array2mat.r")
source("pldag.set.r")
source("ngc.r")
source("grangerLasso.r")
source("grangerTLasso.r")
source("grangerThrLasso.r")
library(glmnet)
library(gglasso)
library(igraph)


###############
## Creates the skeleton of a random sparse d x p x p network 
## with 0.05 sparsity
set.seed(123)
d <- 2
p <- 10
n <- 20
edge = defn_net(d = d, p = p, n = n)
microbenchmark(defn_net(d = d, p = p, n = n), times=10)
# Unit: milliseconds
#    min       lq     mean   median
# 2.45809 2.680687 2.922115 2.831178
#    uq       max     neval
# 3.184873 3.475586    10

###############
## Simulate n iid samples from the granger causality network 
## with p variables observed over T time points
## X is a n x p x T array
T <- 50
error_sd <- 0.2
X <- simulate_data(n, edge, T = T, error_sd = error_sd)
microbenchmark(simulate_data(n, edge, T = T, error_sd = error_sd),
               times=10)
# Unit: milliseconds
#    min       lq    mean   median       uq      max    neval
# 29.51996 31.1318 35.4973 33.23096 36.50567 47.93656    10

###############
## Estimate graphical Granger causality without group effect

# regular
fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.02)
microbenchmark(ngc(X, d = d, method = 'regular', typeIerr = 0.02), 
               times = 10)
# Unit: milliseconds
#   min      lq        mean   median       uq      max    neval
# 39.91263 42.20924 97.60704 44.27759 46.66336 378.1037    10
#

# threshold
fit2 <- ngc(X, d = d, method = 'threshold', typeIerr = 0.02)
microbenchmark(ngc(X, d = d, method = 'threshold', typeIerr = 0.02),
               times = 10)

# Unit: milliseconds
#       min  lq         mean   median     uq      max      neval
#  100.5333 116.0817 125.7991 120.79   127.3765 196.4156    10
# 

# truncate
weights = matrix(1:(d*p)/sum(1:(d*p)), nrow = d, ncol = p)
fit3 <- ngc(X, d = d, method = 'truncate', weights = weights,
            typeIerr = 0.05)
microbenchmark(ngc(X, d = d, method = 'truncate', weights = ,
                   typeIerr = 0.02),
               times = 10)

