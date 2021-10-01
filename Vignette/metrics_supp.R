setwd("/Users/xinyi/Documents/GitHub/ngc/R")
rm(list =ls())
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
library(knitr)

############### Dataset 1
## Creates the skeleton of a random sparse d x p x p network 
## with 0.05 sparsity
set.seed(123)
d <- 2
p <- 16
n <- 20
edge = defn_net(d = d, p = p, n = n, sparsity = 0.1)
grp = c(rep(1,p/2),rep(2,p/2))
edge_grp = defn_net(d = d, p = p, n = n, grp = grp, sparsity = NULL)

phi_true = array(c(edge[2,,],edge[1,,]), dim = c(p,p,d))
phi_true_grp = array(c((edge_grp[2,,]),(edge_grp[1,,])), dim = c(p,p,d))

T <- 50
error_sd <- 0.2


###############
## Perform three different methods with different options for 10 replicates
## dataset has no group effect
sim_num = 10
## Regular lasso
st = Sys.time()
res_lass = list()
for (i in 1:sim_num) {
  set.seed(123*i)
  X <- simulate_data(n, edge, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.05)
  res_lass[[i]] = fit1$estMat
}
end = Sys.time()
t_lass = end - st
t_lass/sim_num

## Threshold lasso
st = Sys.time()
res_tlass = list()
for (i in 1:sim_num) {
  print(i)
  set.seed(123*i)
  X <- simulate_data(n, edge, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'threshold', typeIerr = 0.05)
  res_tlass[[i]] = fit1$estMat
}
end = Sys.time()
t_tlass = end - st
t_tlass/sim_num

## Truncated lasso
st = Sys.time()
res_trlass = list()
for (i in 1:sim_num) {
  print(i)
  set.seed(123*i)
  X <- simulate_data(n, edge, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'truncate', typeIerr = 0.05)
  res_trlass[[i]] = fit1$estMat
}
end = Sys.time()
t_trlass = end - st
t_trlass/sim_num

c(t_lass/sim_num, t_tlass/sim_num, t_trlass/sim_num)
kable(getMeanMedianMetrics(phi_true, res_lass, 'lasso'),"markdown")
kable(getMeanMedianMetrics(phi_true, res_trlass, 'truncated'),"markdown")
kable(getMeanMedianMetrics(phi_true, res_tlass, 'threshold'),"markdown")

###############
## Perform three different methods with different options for 10 replicates
## dataset has group effect
sim_num = 10
## Regular lasso
st = Sys.time()
res_lass = list()
for (i in 1:sim_num) {
  set.seed(123*i)
  X <- simulate_data(n, edge_grp, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.05)
  res_lass[[i]] = fit1$estMat
}
end = Sys.time()
t_lass = end - st
t_lass/sim_num

## Threshold lasso
st = Sys.time()
res_tlass = list()
for (i in 1:sim_num) {
  print(i)
  set.seed(123*i)
  X <- simulate_data(n, edge_grp, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'threshold', typeIerr = 0.05)
  res_tlass[[i]] = fit1$estMat
}
end = Sys.time()
t_tlass = end - st
t_tlass/sim_num

## Truncated lasso
st = Sys.time()
res_trlass = list()
for (i in 1:sim_num) {
  print(i)
  set.seed(123*i)
  X <- simulate_data(n, edge_grp, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'truncate', typeIerr = 0.05)
  res_trlass[[i]] = fit1$estMat
}
end = Sys.time()
t_trlass = end - st
t_trlass/sim_num

c(t_lass/sim_num, t_tlass/sim_num, t_trlass/sim_num)
kable(getMeanMedianMetrics(phi_true_grp, res_lass, 'lasso'),"markdown")
kable(getMeanMedianMetrics(phi_true_grp, res_trlass, 'truncated'),"markdown")
kable(getMeanMedianMetrics(phi_true_grp, res_tlass, 'threshold'),"markdown")

###############
## Perform three different methods with different options for 10 replicates
## dataset has no group effect
sim_num = 10
## Regular lasso
group = c(rep(1,p/2),rep(2,p/2))
st = Sys.time()
res_lass = list()
for (i in 1:sim_num) {
  set.seed(123*i)
  X <- simulate_data(n, edge_grp, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.05, group = group)
  res_lass[[i]] = fit1$estMat
}
end = Sys.time()
t_lass = end - st
t_lass/sim_num

## Threshold lasso
st = Sys.time()
res_tlass = list()
for (i in 1:sim_num) {
  print(i)
  set.seed(123*i)
  X <- simulate_data(n, edge_grp, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'threshold', typeIerr = 0.05, group = group)
  res_tlass[[i]] = fit1$estMat
}
end = Sys.time()
t_tlass = end - st
t_tlass/sim_num


c(t_lass/sim_num, t_tlass/sim_num)
kable(getMeanMedianMetrics(phi_true_grp, res_lass, 'lasso'),"markdown")
kable(getMeanMedianMetrics(phi_true_grp, res_tlass, 'threshold'),"markdown")
