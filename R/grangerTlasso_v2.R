#' Use truncating lasso to estimate graphical Granger causality
#' # note that only work for longitudinal data
#' @param X input array
#' @param d number of time lags to consider
#' @param group group indices
#' @param typeIerr acceptable type I error rate
#' @param typeIIerr acceptable type II error rate
#' @param weights weights for adaptive lasso
#' @param Alasso.power power for adaptive lasso
#' @param eps epsilon used for truncating penalty
#' @param tol tolerance used for BR algorithm
#' @param thresh values smaller than thresh are set to zero
#' @param initialEst initial estimate for BR alg
#' @param maxIter maximum number of iterations for BR alg
#' @return a list including a matrix of estimated coefficients, final lambda value, and time series order estimate
grangerTlasso <-
  function(
    X, 					#input array dim=(n,p,T), last time=Y
    d = NULL, 		#number of time lags to consider
    group = NULL,  
    typeIerr = 0.10, 			#significance level for lambda (...Method=errBased)
    typeIIerr = 0.10,			#acceptable type II error
    weights = weights,  #matrix of weights for Adaptive lasso
    Alasso.power = 1,			#power for Adaptive lasso
    eps = 1e-8,				#epsilon, used for truncating penalty
    tol = 1e-2,				#tolerance, for BR alg
    thresh = 1e-4,		#values smaller than thresh are set to zero
    initialEst=NULL, 	#initial estimate for BR alg
    maxIter=100				#maximum iteration for the BR alg
  ){
    ####### START OF FN #######
    n <- dim(X)[1]
    p <- dim(X)[2]
    tp <- dim(X)[3]

    useAlasso <- !is.null(weights)

    Onep <- matrix(1,p,p)

    newEst <- initialEst
    if (is.null(newEst)){
      newEst <-  array( 0, c(p, p, d) )
    }

    ##scale the X matrix
    for (i in 1:(tp-1)){
      X[,,i] <- scale( X[,,i] )*sqrt(n/(n-1))
    }

    Y <- X[,,tp]		#for compatiblity with previous codes
    YY <- Y

    diff <- tol + 1
    iter <- 0			#number of iterations in the while loop
    FLAG <- TRUE		#flag for whether to stop the while loop
    d <- d		#number of lags to consider in estimation
    lambda <- NULL
    
    #while loop for BR algorithm
    #cat("BR alg started! \n")
    while ( (diff > tol) && (iter<=maxIter) && (FLAG==TRUE) ){
      print(diff)
      print(FLAG)
      oldEst <- newEst
      iter <- iter + 1
      #cat("BR loop No: ", iter, "\n")

      ##estimation loop
      jj <- 1
      CONTINUE <- TRUE		#flag for continuation of inner estimation loop
      while (jj <= d){
        print(jj)
        XX <- X[,,(tp-jj)]
        
        # #Find the residual as Y # ?(similar to Gauss-Sidel)
        # #? Error in array2mat(X[, , -c(tp - jj, tp)]) %*% theta[i, ] :
        # #? non-conformable arguments
        # theta <- array2mat(newEst[,,-(tp-jj)])
        # for (i in 1:p){
        #   YY[,i] <- Y[,i] -
        #     array2mat(X[, , -c(tp-jj,tp)]) %*% theta[i,]
        # }
        
        #? Change to below
        # ddddd = array(1:12, dim = c(2,2,3))
        # use this matrix(ddddd, nrow = dim(ddddd)[1], byrow = FALSE)
        # instead of array2mat
        theta <- array2mat(newEst[,,-(jj)])
        lag_ind = seq(tp-1, tp-d, -1)
        for (i in 1:p){
          YY[,i] <- Y[,i] -
            array2mat(X[, , rev(lag_ind[-jj])]) %*% theta[i,]
        }

        # ## calculate truncating factor
        # if (jj>1){
        #   ##? previousNNZ <- sum( abs(newEst[,,(tp-jj+1)]) > thresh )/(p^2)
        #   previousNNZ <- sum( abs(newEst[,,(d-jj+2)]) > thresh )/(p^2)
        #   if( previousNNZ < (typeIIerr/(tp-jj)) ){
        # 
        #     newEst[,,1:(d-jj+1)] <- matrix(0,p,p*(d-jj+1))
        #     # d <- jj - 1
        #     CONTINUE <- FALSE
        #   }
        # }

        ## calculate est
        if (CONTINUE==TRUE){
          if (!useAlasso){
            tempp <- pldag.set(X1=XX, X2=YY, group = group,
                               sigLevel=typeIerr, useWghts=FALSE,
                               wantScale=TRUE, useTLASSO=TRUE, d=d)

            ##? newEst[,,(tp-jj)] <- tempp$AA
            newEst[,, (d-jj+1)] <- tempp$AA
          }else{
            print(dim(weights))
            W <- (abs(weights[,,(d-jj+1)]) + eps*Onep)^(-Alasso.power)
            if (sum(W < 1) > 0){
              W[(W < 1)] <- 1
            }
            tempp <-pldag.set(XX, YY, group = group,
                              sigLevel=typeIerr, useWghts=TRUE,
                              wghts=W, wantScale=TRUE, useTLASSO=TRUE, d=d)

            newEst[,,(d-jj+1)] <- tempp$AA
          }#endif(!useAlasso)
        }#endif(CONTINUE)
        jj <- jj + 1

      }#endwhile(jj)

      diff <- sum( abs(newEst-oldEst) ) / sum( abs(newEst)>eps )
      if ((diff == Inf) || is.na(diff))
      {
        diff = tol + 1
      }
    }#endwhile()
    lambda <- tempp$lambda
    intercepts <- tempp$intercepts

    if (diff > tol){
      cat("WARNING: BR Alg not converged! Increase maxIter. \n")
    }

    return(list(estMat=newEst, lambda=lambda, tsOrder=jj-1, intercepts = intercepts))
  }

