setwd("/Users/xinyi/Documents/GitHub/network_ngc/R")
rm(list =ls())
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
n <- 30
edge = defn_net(d = d, p = p, n = n, sparsity = NULL)
grp = c(rep(1,5),rep(2,5))
set.seed(123)
edge_grp = defn_net(d = d, p = p, n = n, grp = grp, sparsity = NULL)
# microbenchmark(defn_net(d = d, p = p, n = n), times=10)

phi_true = array(c(edge[2,,],edge[1,,]), dim = c(p,p,d))
phi_true_grp = array(c((edge_grp[2,,]),(edge_grp[1,,])), dim = c(p,p,d))

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
# microbenchmark(simulate_data(n, edge, T = T, error_sd = error_sd),
               # times=10)

# Unit: milliseconds
#    min       lq    mean   median       uq      max    neval
# 29.51996 31.1318 35.4973 33.23096 36.50567 47.93656    10

###############
## Estimate graphical Granger causality without group effect
# regular
fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.05)
# set.seed(123)
# microbenchmark(ngc(X, d = d, method = 'regular', typeIerr = 0.02), 
#                times = 10)
xtable(checkEstimation(phi_true, fit1$estMat, method_name = "Lasso")$metrics,)

fit1 = ngc(X, d = d, method = 'regular')
xtable(checkEstimation(phi_true, fit1$estMat, method_name = "Lasso")$metrics,)

# Unit: milliseconds
#   min      lq        mean   median       uq      max    neval
# 39.91263 42.20924 97.60704 44.27759 46.66336 378.1037    10
#

# regular with grouping by matrices
group = c(rep(1,5),rep(2,5))
fit12 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = group)
# microbenchmark(ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = group), 
#                times = 10)

# Unit: milliseconds
#   min      lq        mean   median       uq      max    neval
#  86.7227 87.88322 97.08153 93.91508 105.145 113.3507    10
xtable(checkEstimation(phi_true, fit12$estMat, method_name = "Lasso")$metrics)

# regular with grouping by time
fit13 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = group, groupByTime = TRUE)
xtable(checkEstimation(phi_true, fit13$estMat, method_name = "Lasso")$metrics)

# regular with grouping by time without grouping by states
fit13 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, groupByTime = TRUE)
# microbenchmark(ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = group, groupByTime = TRUE), 
               # times = 10)
# Unit: milliseconds
#   min      lq        mean   median       uq      max    neval
#  112.6255 118.404 125.5603 122.6254 130.919  152.1969    10
xtable(checkEstimation(phi_true, fit13$estMat, method_name = "Lasso")$metrics)


# threshold
fit2 <- ngc(X, d = d, method = 'threshold', typeIerr = 0.02)
microbenchmark(ngc(X, d = d, method = 'threshold', typeIerr = 0.02),
               times = 10)
dim(fit2$estMat)
xtable(checkEstimation(phi_true, fit2$estMat, method_name = "Lasso")$metrics)

# Unit: milliseconds
#       min  lq         mean   median     uq      max      neval
#  100.5333 116.0817 125.7991 120.79   127.3765 196.4156    10
# 

# truncate
d = 5
# weights = matrix(1:(d*p)/sum(1:(d*p)), nrow = d, ncol = p)
# fit3 <- ngc(X, d = d, method = 'truncate', weights = weights,
#             typeIerr = 0.05)
d = 2
fit3 <- ngc(X, d= 2,method = 'truncate', 
            typeIerr = 0.05)
microbenchmark(ngc(X, method = 'truncate', 
                   typeIerr = 0.05),
               times = 10)
xtable(checkEstimation(phi_true, fit3$estMat[,,dim(fit3$estMat)[3]-(1:2)], method_name = "Lasso")$metrics)

# Unit: milliseconds
# min       lq      mean      median    uq        max       neval
# 915.2179 1037.183 1065.971 1053.627  1080.478 1342.806    10
# 

microbenchmark(ngc(X, d = 5, method = 'truncate', 
                  typeIerr = 0.05),
               times = 10)
# Unit: milliseconds
# min       lq       mean      median      uq      max        neval
# 454.4371 472.3171 480.4898  480.6106  497.928    504.0195    10

### functions for simulation

ngc_sim1 = function(n=20, d=2, p=10, T=50, group = NULL, groupByTime = FALSE, seedj) {
  set.seed(seedj)
  edge = defn_net(d = d, p = p, n = n, sparsity = 0.1)
  error_sd <- 0.2
  X <- simulate_data(n, edge, T = T, error_sd = error_sd)
  fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = group, groupByTime = groupByTime)
  return(list(edge, fit1))
}
ngc_sim2 = function(n=20, d=2, p=10, T=50, group = NULL, groupByTime = FALSE, seed) {
  set.seed(seed)
  edge = defn_net(d = d, p = p, n = n, sparsity = 0.1)
  error_sd <- 0.2
  X <- simulate_data(n, edge, T = T, error_sd = error_sd)
  fit2 <- ngc(X, d = d, method = 'threshold', typeIerr = 0.02, group = group, groupByTime = groupByTime)
  return(list(edge, fit2))
}

ngc_sim3 = function(n=20, d=2, p=10, T=50, group = NULL, groupByTime = FALSE, seed) {
  set.seed(seed)
  edge = defn_net(d = d, p = p, n = n, sparsity = 0.1)
  error_sd <- 0.2
  X <- simulate_data(n, edge, T = T, error_sd = error_sd)
  fit3 <- ngc(X, d=d, method = 'truncate', 
              typeIerr = 0.02)
  return(list(edge, fit3))
}

ngc_sim4 = function(n=20, d=2, p=10, T=50, group = NULL, groupByTime = FALSE, seed) {
  set.seed(seed)
  edge = defn_net(d = d, p = p, n = n, sparsity = 0.1)
  error_sd <- 0.2
  X <- simulate_data(n, edge, T = T, error_sd = error_sd)
  fit3 <- ngc(X, method = 'truncate', 
              typeIerr = 0.02)
  return(list(edge, fit3))
}

## Check n
sim_num = 10
n_set = c(20, 50, 100, 200)
time_nall = matrix(0, nrow = 4, ncol = length(n_set))
group = rep(1:2,each = 5)

for (i in 1:length(n_set)) {
  print(i)
  time_ls = matrix(0, nrow = sim_num, ncol = 4)
  error_ls = rep(list(matrix(0, nrow = sim_num, ncol = 4)), 4)
  for (j in 1:sim_num) {
    set.seed(j*12345)
    p = 10
    d = 2
    n = n_set[i]
    edge = defn_net(d = d, p = p, n = n, sparsity = 0.1)
    error_sd <- 0.2
    X <- simulate_data(n, edge, T = T, error_sd = error_sd)
    phi_true = array(c(edge[2,,],edge[1,,]), dim = c(p,p,d))
    
    print(paste("J:",j))
    st = Sys.time()
    tmp1 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = NULL, groupByTime = FALSE)
    end = Sys.time()
    time_ls[j, 1] = end-st
    error_ls[[1]][j,] = checkEstimation(phi_true, tmp1$estMat, method_name = "Lasso")$metrics[2,]
    
    st = Sys.time()
    tmp12 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = rep(1:2,each = 5), groupByTime = FALSE)
    end = Sys.time()
    time_ls[j, 2] = end-st
    error_ls[[2]][j,] = checkEstimation(phi_true, tmp12$estMat, method_name = "Lasso")$metrics[2,]
    
    st = Sys.time()
    tmp13 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = NULL, groupByTime = TRUE)
    end = Sys.time()
    time_ls[j, 3] = end-st
    error_ls[[3]][j,] = checkEstimation(phi_true, tmp13$estMat, method_name = "Lasso")$metrics[2,]
    
    st = Sys.time()
    tmp14 = ngc(X, d = d, method = 'regular', typeIerr = 0.02, group = rep(1:2,each = 5), groupByTime = TRUE)
    end = Sys.time()
    time_ls[j, 4] = end-st
    error_ls[[4]][j,] = checkEstimation(phi_true, tmp14$estMat, method_name = "Lasso")$metrics[2,]
    
  }
  
  
  time_nall[, i] = apply(time_ls,2,median)
}

n_set = c(20, 50, 100, 200)
time_nall = matrix(0, nrow = 4, ncol = length(n_set))
for (i in 1:length(n_set)) {
  tmp1 = microbenchmark(ngc_sim1(n = n_set[i]), times = 10)
  tmp12 = microbenchmark(ngc_sim1(n = n_set[i]), times = 10, group = group)
  tmp13 = microbenchmark(ngc_sim1(n = n_set[i], group = group, groupByTime = TRUE),  times = 10)
  d = 2
  weights = matrix(1:(d*p)/sum(1:(d*p)), nrow = d, ncol = p)
  tmp14 = microbenchmark(ngc_sim1(n = n_set[i], group = group, groupByTime = TRUE, weights = weights), times = 10)
  tmp2 = microbenchmark(ngc_sim2(n = n_set[i]), times = 10)
  tmp3 = microbenchmark(ngc_sim3(n = n_set[i]), times = 10)
  tmp4 = microbenchmark(ngc_sim4(n = n_set[i]), times = 10)
  time_nall[1, i] = mean(tmp1[,2])/1e9
  time_nall[2, i] = mean(tmp2[,2])/1e9
  time_nall[3, i] = mean(tmp3[,2])/1e9
  time_nall[4, i] = mean(tmp4[,2])/1e9
}

n_set = c(20, 50, 100, 200)
time_nall = matrix(0, nrow = 4, ncol = length(n_set))
for (i in 1:length(n_set)) {
  tmp1 = microbenchmark(ngc_sim1(n = n_set[i]), times = 10)
  tmp12 = microbenchmark(ngc_sim1(n = n_set[i]), times = 10, group = group)
  tmp13 = microbenchmark(ngc_sim1(n = n_set[i], group = group, groupByTime = TRUE),  times = 10)
  d = 2
  weights = matrix(1:(d*p)/sum(1:(d*p)), nrow = d, ncol = p)
  tmp14 = microbenchmark(ngc_sim1(n = n_set[i], group = group, groupByTime = TRUE, weights = weights), times = 10)
  tmp2 = microbenchmark(ngc_sim2(n = n_set[i]), times = 10)
  tmp3 = microbenchmark(ngc_sim3(n = n_set[i]), times = 10)
  tmp4 = microbenchmark(ngc_sim4(n = n_set[i]), times = 10)
  time_nall[1, i] = mean(tmp1[,2])/1e9
  time_nall[2, i] = mean(tmp2[,2])/1e9
  time_nall[3, i] = mean(tmp3[,2])/1e9
  time_nall[4, i] = mean(tmp4[,2])/1e9
}

n_set = c(20, 50, 100, 200)
time_nall = matrix(0, nrow = 4, ncol = length(n_set))
for (i in 1:length(n_set)) {
  tmp1 = microbenchmark(ngc_sim1(n = n_set[i]), times = 10)
  tmp12 = microbenchmark(ngc_sim1(n = n_set[i]), times = 10, group = group)
  tmp13 = microbenchmark(ngc_sim1(n = n_set[i], group = group, groupByTime = TRUE),  times = 10)
  d = 2
  weights = matrix(1:(d*p)/sum(1:(d*p)), nrow = d, ncol = p)
  tmp14 = microbenchmark(ngc_sim1(n = n_set[i], group = group, groupByTime = TRUE, weights = weights), times = 10)
  tmp2 = microbenchmark(ngc_sim2(n = n_set[i]), times = 10)
  tmp3 = microbenchmark(ngc_sim3(n = n_set[i]), times = 10)
  tmp4 = microbenchmark(ngc_sim4(n = n_set[i]), times = 10)
  time_nall[1, i] = mean(tmp1[,2])/1e9
  time_nall[2, i] = mean(tmp2[,2])/1e9
  time_nall[3, i] = mean(tmp3[,2])/1e9
  time_nall[4, i] = mean(tmp4[,2])/1e9
}

## Check p
p_set = c(5, 10, 20, 50)
time_pall = matrix(0, nrow = 4, ncol = length(p_set))

for (i in 1:length(p_set)) {
  tmp1 = microbenchmark(ngc_sim1(p = p_set[i]), times = 10)
  tmp2 = microbenchmark(ngc_sim2(p = p_set[i]), times = 10)
  tmp3 = microbenchmark(ngc_sim3(p = p_set[i]), times = 10)
  tmp4 = microbenchmark(ngc_sim4(p = p_set[i]), times = 10)
  time_pall[1, i] = median(tmp1[,2])/1e9
  time_pall[2, i] = median(tmp2[,2])/1e9
  time_pall[3, i] = median(tmp3[,2])/1e9
  time_pall[4, i] = median(tmp4[,2])/1e9
}

## Check T
T_set = c(50, 100, 150, 200)
time_Tall = matrix(0, nrow = 4, ncol = length(T_set))

for (i in 1:length(T_set)) {
  tmp1 = microbenchmark(ngc_sim1(T = T_set[i]), times = 10)
  tmp2 = microbenchmark(ngc_sim2(T = T_set[i]), times = 10)
  tmp3 = microbenchmark(ngc_sim3(T = T_set[i]), times = 10)
  tmp4 = microbenchmark(ngc_sim4(T = T_set[i]), times = 10)
  time_Tall[1, i] = mean(tmp1[,2])/1e9
  time_Tall[2, i] = mean(tmp2[,2])/1e9
  time_Tall[3, i] = mean(tmp3[,2])/1e9
  time_Tall[4, i] = mean(tmp4[,2])/1e9
}

## Check d



## Plot
par(mfrow = c(1,3))
plot(n_set, time_nall[1,], ylim = c(min(time_nall), max(time_nall)),
     type = "l", main = "Time Change with n (mean)", ylab = "Time (Second)", xlab = "n")
for (i in 2:4){
  lines(n_set, time_nall[i,], col = i)
}

plot(p_set, time_pall[1,], ylim = c(min(time_pall), max(time_pall)),
     type = "l", main = "Time Change with p (mean)", ylab = "Time (Second)", xlab = "p")
for (i in 2:4){
  lines(p_set, time_pall[i,], col = i)
}

plot(T_set, time_Tall[1,], ylim = c(min(time_Tall), max(time_Tall)),
     type = "l", main = "Time Change with T (mean)", ylab = "Time (Second)", xlab = "T")
for (i in 2:4){
  lines(T_set, time_Tall[i,], col = i)
}


## estimate metrics
sim_num = 10
set.seed(123)
phi_est_all = list()
for (i in 1:sim_num) {
  fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.02)
  phi_est_all[[i]] = fit1$estMat
}
xtable(getMeanMedianMetrics(phi_true, phi_est_all) )

               