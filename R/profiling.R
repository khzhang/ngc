setwd("/Users/xinyi/Documents/GitHub/network_ngc/R")
library('microbenchmark')
source("defn.R")
source("simulate.r")
source("array2mat.r")
source("pldag.set.r")
source("ngc.r")
source("grangerLasso.r")
source("grangerTLasso_v3.r")
source("grangerThrLasso.r")
source("metrics.R")
library(glmnet)
library(gglasso)
library(igraph)
library(xtable)


###############
## Creates the skeleton of a random sparse d x p x p network 
## with 0.05 sparsity
set.seed(123)
d <- 2
p <- 10
n <- 20
edge = defn_net(d = d, p = p, n = n, sparsity = 0.1)
set.seed(123)
grp = c(rep(1,5),rep(2,5))
edge_grp = defn_net(d = d, p = p, n = n, grp = grp, sparsity = 0.1)
microbenchmark(defn_net(d = d, p = p, n = n), times=10)
phi_true = array(c((edge[1,,]),(edge[2,,])), dim = c(10,10,2))
phi_true_grp = array(c((edge_grp[1,,]),(edge_grp[2,,])), dim = c(10,10,2))

# Unit: milliseconds
#    min       lq     mean   median
# 2.45809 2.680687 2.922115 2.831178
#    uq       max     neval
# 3.184873 3.475586    10

###############
## Simulate n iid samples from the granger causality network 
## with p variables observed over T time points
## X is a n x p x T array
set.seed(123)
T <- 300
error_sd <- 0.2
X <- simulate_data(n, edge, T = T, error_sd = error_sd)
X_grp <- simulate_data(n, edge_grp, T = T, error_sd = error_sd)

Rprof("reg.out")  
fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.05, group = group, groupByTime = TRUE)
plot.ngc(fit1, ngc.type = "dag")
fit1_pred = predict(fit1, 2)
Rprof(NULL)
summaryRprof("reg.out")


Rprof(filename = "reg_cv.out")  
fit1 = ngc(X, d=d, method = 'regular', group = group, groupByTime = TRUE)
plot.ngc(fit1, ngc.type = "dag")
fit1_pred = predict(fit1, 2)
Rprof(NULL)
tmp = summaryRprof("reg_cv.out")


library(profvis)
group = c(rep(1,5),rep(2,5))
profvis({
  fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.05, group = group, groupByTime = TRUE)
  plot.ngc(fit1, ngc.type = "dag")
  fit1_pred = predict(fit1, 2)
})

profvis({
  fit1 = ngc(X, d=d, method = 'regular', group = group, groupByTime = TRUE)
  plot.ngc(fit1, ngc.type = "dag")
  fit1_pred = predict(fit1, 2)
})

group = c(rep(1,5),rep(2,5))
profvis({
  fit1 = ngc(X, d = d, method = 'threshold', typeIerr = 0.05, group = group, groupByTime = TRUE)
  plot.ngc(fit1, ngc.type = "dag")
  fit1_pred = predict(fit1, 2)
})

profvis({
  fit1 = ngc(X, d=d, method = 'threshold', group = group, groupByTime = TRUE)
  plot.ngc(fit1, ngc.type = "dag")
  fit1_pred = predict(fit1, 2)
})

