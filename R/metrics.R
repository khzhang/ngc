
# phi_est = fit1$estMat
# phi_true = array(c(edge[1,,],edge[2,,]), dim = c(10,10,2))

checkEstimation = function(phi_true, phi_est, method_name = "Lasso" ){
  ## True edges
  N = length(phi_true) 
  count_nonzero = sum(unlist(phi_true) != 0)
  count_zero = N - count_nonzero
  
  est_error =  matrix(0, 3, 7)
  est_error[3,1] = method_name
  est_error[1,2] = "REE"
  est_error[1,3] = "Fobinius Norm" 
  est_error[1,4] = "TZ"
  est_error[2,4] = count_zero
  est_error[2,1] = "TRUTH"
  est_error[1,5] = "TNZ"
  est_error[2,5] = count_nonzero
  est_error[1,6] = "FZ"
  est_error[1,7] = "FNZ"

  
  ## Estimation
  count_falsezero = sum(c(phi_true != 0 & phi_est == 0))
  count_falsenonzero = sum(c(phi_true == 0 & phi_est != 0))
  count_truezero = sum(c(phi_true == 0 & phi_est == 0))
  count_truenonzero = sum(c(phi_true != 0 & phi_est != 0))
  l2_norm = sqrt(sum((phi_true - phi_est)^2))
  REE = l2_norm/sum(phi_true^2)
  
  ## Metrics
  est_error[3,2] = REE 
  est_error[3,3] = l2_norm 
  est_error[3,4] = count_truezero
  est_error[3,5] = count_truenonzero
  est_error[3,6] = count_falsezero
  est_error[3,7] = count_falsenonzero
  
  est_error[2:3,2:7] = round(as.numeric(est_error[2:3,2:7]), digits = 3)
  FPR = as.numeric(est_error[3,7])/(as.numeric(est_error[3,7])+as.numeric(est_error[3,4]))
  TPR = as.numeric(est_error[3,5])/(as.numeric(est_error[3,5])+as.numeric(est_error[3,6]))
  est_error2 = est_error[c(1,3),c(2,3)]
  est_error2 = cbind(est_error2,c("FPR",round(FPR, 3)),c("TPR",round(TPR, 3)))
  
  return( list(REE = REE, 
               l2_norm = l2_norm,
               true_zero = count_truezero,
               true_nonzero = count_truenonzero, 
               false_zero = count_falsezero,
               false_nonzero = count_falsenonzero, 
               metrics = est_error2) )
}

getMeanMedianMetrics = function(phi_true, phi_est_all) {
  ## True edges
  N = length(phi_true) 
  count_nonzero = sum(unlist(phi_true) != 0)
  count_zero = N - count_nonzero
  
  est_error =  matrix(0, 3, 9)
  est_error[3,1] = method_name
  est_error[1,2] = "REE"
  est_error[1,3] = "SD(REE)" 
  est_error[1,4] = "Fobinius Norm"
  est_error[1,5] = "SD(Fobinius Norm)" 
  
  est_error[1,6] = "TZ"
  est_error[2,6] = count_zero
  est_error[2,1] = "TRUTH"
  est_error[1,7] = "TNZ"
  est_error[2,7] = count_nonzero
  est_error[1,8] = "FZ"
  est_error[1,9] = "FNZ"
  
  ## Estimation
  N = length(phi_est_all) 
  l2_norm = rep(0,N) 
  REE = rep(0,N)
  true_zero = rep(0,N) 
  true_nonzero = rep(0,N)
  false_zero = rep(0,N)
  false_nonzero = rep(0,N)
  for (i in 1:N) {
    phi_temp = phi_est_all[[i]]
    res = checkEstimation(phi_true, phi_temp)
    l2_norm[i] = res$l2_norm
    REE[i] = res$REE
    true_zero[i] = res$true_zero
    true_nonzero[i] = res$true_nonzero
    false_zero[i] = res$false_zero
    false_nonzero[i] = res$false_nonzero
  }
  
  REE_mean = mean(REE)
  REE_sd = sd(REE)
  l2_norm_mean = mean(l2_norm)
  l2_norm_sd = sd(l2_norm)
  true_zero_median = median(true_zero)
  true_nonzero_median = median(true_nonzero)
  false_zero_median = median(false_zero)
  false_nonzero_median = median(false_nonzero)
  
  ## Metrics
  est_error[3,2] = REE_mean 
  est_error[3,3] = REE_sd 
  est_error[3,4] = l2_norm_mean 
  est_error[3,5] = l2_norm_sd 
  
  est_error[3,6] = count_truezero
  est_error[3,7] = count_truenonzero
  est_error[3,8] = count_falsezero
  est_error[3,9] = count_falsenonzero
  
  est_error[2:3,2:9] = round(as.numeric(est_error[2:3,2:9]),digits = 3)
  FPR = count_falsenonzero/(count_falsenonzero+count_truezero)
  TPR = count_truenonzero/(count_truenonzero+count_falsezero)
  est_error2 = est_error[c(1,3),2:5]
  est_error2 = cbind(est_error2, c("FPR",round(FPR,3)), c("TPR",round(TPR,3)))
  
  return(est_error2)
}

# xtable(getMeanMedianMetrics)
