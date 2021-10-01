## install packages
install.packages("devtools")
devtools::install_github("khzhang/ngc", build_vignettes=T)

library(Matrix)
library(glmnet)
library(gglasso)
library(igraph)
library(xtable)

## Simulate data
###############
## Creates the skeleton of a random sparse d x p x p network 
## with 0.1 sparsity

set.seed(123)
d <- 2
p <- 16
n <- 10
edge <- defn_net(d = d, p = p, n = n, sparsity = 0.1)

## Simulate n iid samples from the granger causality network 
## with p variables observed over T time points
## X is a n x p x T array
T <- 100
error_sd <- 0.2
X <- simulate_data(n, edge, T = T, error_sd = error_sd)

## Estimations and Visualization
###############
## Estimate graphical Granger causality using regular lasso approach
fit1 = ngc(X, d = d, method = 'regular', typeIerr = 0.05)
fit1_typeI = ngc(X, d=2, typeIerr = 0.02)
fit1_typeI$estMat
fit1_typeI$lambda

## Plot the Granger causality
plot.ngc(fit1, ngc.type = "dag")

## Predict the next time point
predict(fit1, 1)


