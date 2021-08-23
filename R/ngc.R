#' Estimate graphical Granger causality
#' @param X input array
#' @param d number of time lags to consider
#' @param method estimation method to use; options are "regular", "truncate", or "threshold"
#' @param group vector of group indices of length p or p*d; if null, no group structure
#' @param groupByTime whether to group covariates across time points
#' @param typeIerr acceptable type I error rate
#' @param typeIIerr acceptable type II error rate
#' @param weights weights for adaptive lasso
#' @param thresholdConstant constant used to compute threshold
#' @param refit whether to refit a linear regression after initial thresholding
#' @param edgeThreshold absolute value threshold for including edge in graph
#' @param covNames covariate names
#' @return a list including a matrix of estimated coefficients, final lambda value, and time series order estimate
#' @export ngc

ngc <-
  function(
    X, #input array dim=(n,p,T) (longitudinal), or (p,T) (time series); last time=Y %%what is Y?
    d = NULL, #number of time lags to  regularconsider
    method = 'regular', #method to use. options are "regular", "truncate", and "threshold"
    group = NULL, #vector of group indices of length p or p*d; if null, no group structure
    groupByTime = FALSE, #whether to group covariates over time 
    typeIerr = NULL, #acceptable type I error rate; if provided, error-based lasso is fitted
    typeIIerr = 0.1, #acceptable type II error rate
    weights = NULL, #matrix of weights for adaptive lasso. If no weights are provided, use regular lasso.
    thresholdConstant = NULL, #constant used for calculating threshold value
    refit = FALSE, #whether to refit a linear regression after initial thresholding
    covNames = NULL #covariate names
  ){
    ####### START OF FN #######
    if (method != 'regular' & method != 'truncate' & method != 'threshold')
    {
      stop('Invalid estimation method specified')
    }

    if (!is.array(X) | length(dim(X)) <2 | length(dim(X)) > 3)
    {
      stop('Invalid X')
    }

    if (length(dim(X))==2)
    { 
      # size of time series data
      n <- 1
      p <- dim(X)[1]
      len <- dim(X)[2]
    }
    else
    { 
      # size of longitudinal data
      n <- dim(X)[1]
      p <- dim(X)[2]
      len <- dim(X)[3]
    }

    if (is.null(d))
    {
      # use the maximum possible lags
      d <- len-1
    }

    if (d >= len)
    {
      stop('Number of time lags to consider cannot exceed number of time points')
    }

    #Set up replicates for the time series case
    #Put X into array format
    #The transformed matrix has (len-d) replicates over d+1 time points
    if (n == 1)
    {
      if (d >= len-1)
      {
        stop('Number of time lags to consider must be restricted in order to fit time series')
      }
      cat('Warning: stationarity assumption is required for time series data')
      xMat <- X
      n <- len-d
      len <- d+1
      X <- array(0, c(n,p,len))
      for (i in 1:n)
      {
        X[i,,] <- xMat[,i:(d+i)] #%% include current time point
      }
    }
    
    if (!is.null(covNames))
    {
      if (length(covNames) != p)
      {
        stop("Number of covariate names must match number of covariates")
      }
    }
    
    if (!is.null(group))
    {
      if (length(group)!=p)
      {
        stop('Invalid group specification')
      }
      if (!is.numeric(group))
      {
        stop('Groups must be specified with consecutive integers')
      }
      if (!all.equal(order(group), 1:p))
      {
        stop("Groups must be specified with consecutive integers")
      }
      # apply groups across time points 
      if (groupByTime)
      {
        group <- rep(group, d)
      }
      else
      {
        ngrp = length(unique(group))
        group <- group + rep(seq(0, (d-1)*ngrp, by = ngrp), each = p)
      }
    } 
    
    if (groupByTime)
    {
      group <- rep(1:p, d)
    }

    if (method == 'regular')
    {
      fit <- grangerLasso(X, d = d, group = group, typeIerr = typeIerr,
                          weights = weights)
    }

    else if (method == 'truncate')
    {
      fit <- grangerTlasso(X, d = d, group = group, typeIerr = typeIerr,
                           typeIIerr = typeIIerr, weights = weights)
    }

    else #threshold
    {
      fit <- grangerThrLasso(X, d = d, group = group, typeIerr = typeIerr,
                             typeIIerr = typeIIerr, weights = weights,
                             thresholdConstant = thresholdConstant,
                             refit = refit)
    }
    
    
    edgeIx <- which(fit$estMat != 0, arr.ind = T)
    edgeCount <- dim(edgeIx)[1]
    if (is.null(fit$tsOrder))
    {
      tsOrder <- ifelse(edgeCount > 0, max(d-edgeIx[,3]+1), 0)
      fit$tsOrder <- ifelse(!is.null(tsOrder), tsOrder, 0)
    }
    
    if (fit$tsOrder == 0){
      dagMat <- Matrix(0, nrow=p*(d+1), ncol=p*(d+1), sparse = TRUE)
      ringMat <- Matrix(0, nrow=p, ncol=p)
      edgeIx2 = edgeIx
    } else {
      edgeIx2 <- (which(fit$estMat[,,dim(fit$estMat)[3]+1-seq(fit$tsOrder,1,-1)] != 0, arr.ind = T))
      if (ncol(edgeIx2) != 3) {
        edgeIx2 <- cbind(edgeIx2, 1)
      }
      dagMat <- Matrix(0, nrow=p*(fit$tsOrder+1), ncol=p*(fit$tsOrder+1), sparse = TRUE)
      ringMat <- Matrix(0, nrow=p, ncol=p)
    }
    
    # dagMat <- Matrix(0, nrow=p*(d+1), ncol=p*(d+1), sparse = TRUE)
    # ringMat <- Matrix(0, nrow=p, ncol=p)
    
    if (edgeCount > 0)
    {
      for (i in 1:edgeCount)
      {
        # print(i)
        edge <- edgeIx[i,]
        edge2 <- edgeIx2[i,]
        pStart <- edge[2]
        pEnd <- edge[1]
        lag <- edge[3]
        
        dagMat[((edge2[3]-1)*p + pStart),(fit$tsOrder*p + pEnd)] <- fit$estMat[pEnd, pStart, lag]
        ringMat[pStart, pEnd] <- ringMat[pStart, pEnd] + fit$estMat[pEnd, pStart, lag]^2
      } 
    }
    
    fit$dag <- graph_from_adjacency_matrix(dagMat, mode = 'directed', weighted = TRUE)
    fit$ring <- graph_from_adjacency_matrix(sqrt(ringMat), mode = 'directed', weighted = TRUE)
    fit$method <- method
    fit$n <- n
    fit$p <- p
    fit$len <- len
    fit$d <- d
    fit$group <- group
    fit$X <- X
    fit$covNames <- covNames
    class(fit) <- "ngc"
    return(fit)
  }

#' Plot DAG of network or ring graph showing Granger causality
#' @param fit object of class ngc
#' @param ngc.ring whether to plot ring graph
plot.ngc <- 
  function(
    fit, #object of class ngc
    ngc.type = "dag", #"dag" or "granger"
  ...){
    if (class(fit) != "ngc")
    {
      stop("Class of argument must be ngc")
    }
    p <- fit$p
    d <- fit$tsOrder
    covNames <- fit$covNames
    group <- fit$group
    if (ngc.type == "granger")
    {
      g <- fit$ring
      if (is.null(E(g)$weight))
      {
        edgeThickness = 0
      }
      else
      {
        edgeThickness = scale(E(g)$weight^2/mean(E(g)$weight^2),scale = T)*100
      }
      edgeThickness <- ifelse(edgeThickness > 0.5, edgeThickness, 0.5)
      edgeThickness <- ifelse(edgeThickness < 2, edgeThickness, 2)
      labelCex <- max(min(10/p, 1), 0.3)
      aRatio <- (d/p)/2
      arrowSize <- 0.3*labelCex
      par(mar=c(0.5, 0.5, 0.5, 1.5))
      plot(g, layout = layout_in_circle(g), asp = aRatio,edge.arrow.size = arrowSize,
           vertex.shape = "none", cex = 1.2, edge.width = edgeThickness, ...)
    }
    else
    {
      xcoords = rep(1:(d+1), each=p)
      ycoords = rep(p:1, d+1)
      g <- fit$dag
      if (d != 0 ){
        layout_matrix = matrix(c(xcoords, ycoords), ncol=2)
        
        groupList <- NULL
        if (!is.null(group))
        {
          groupList <- lapply(unique(group),function(x){which(group==x)})
        }
        par(mar=c(1.5, 0.5, 0.5, 1.5))
        edgeColor = ifelse(E(g)$weight > 0, "blue", "red")
        if (is.null(E(g)$weight))
        {
          edgeThickness = 0
        }
        else
        {
          edgeThickness <- scale(E(g)$weight^2/mean(E(g)$weight^2),scale = T)*100
        }
        #control maximum and minimum thickness
        edgeThickness <- ifelse(edgeThickness > 0.5, edgeThickness, 0.5)
        edgeThickness <- ifelse(edgeThickness < 2, edgeThickness, 2)
        labelCex <- max(min(10/p, 1), 0.3)
        arrowSize <- 0.3*labelCex
        #curve edges that are more than 1 lag
        edgeTails <- tail_of(g, E(g))
        edgeCurvature <- (edgeTails <= p*(d-1))*0.25
        edgeCurvature <- edgeCurvature*(-1)^((head_of(g, E(g)) %% p) < (edgeTails %% p))
        aRatio <- (d/p)/2
        plot(g, asp = aRatio, layout = layout_matrix,
             mark.groups = groupList, mark.border = NA,
             vertex.label.cex = labelCex * 0.6,
             vertex.label = rep(1:p, d+1), vertex.shape = "none",
             edge.color = edgeColor, edge.width = edgeThickness,
             edge.arrow.size = arrowSize, edge.curved = edgeCurvature,
             rescale = FALSE, xlim = c(1, d+1), ylim = c(0, p), ...)
        text(0, -0.5, "Lag", cex = labelCex)
        lagStep <- ifelse(d < 10, 1, 5)
        for (i in seq(lagStep, d, lagStep))
        {
          text(i, -0.5, d-i+1, cex = labelCex*0.6)
        }
        if (!is.null(covNames))
        {
          legend(d+1+2*aRatio, p+0.5, paste(1:p, covNames, sep = " - "), cex = labelCex, ncol = p%/%10 + 1, title = "Legend")
        }
      } else {
        plot(g,...)
      }
    }
  }

#' Predict covariate values at a given time point
#' @param fit object of class ngc from which to predict
#' @param tp time point at which to predict covariate values 
#' e.g. if tp = 2, output is fitted covariates at time T+2
#' @return n x p matrix of fitted covariates for each replicate
predict.ngc <- 
  function(
    fit, #object of class ngc
    tp = 0#time point at which to predict covariate value
  ){
    if (class(fit) != "ngc")
    {
      stop("Class of argument must be ngc")
    }
    if (tp < 0)
    {
      stop("Specify current or future time point")
    }
    n <- fit$n
    p <- fit$p
    d <- fit$d
    len <- fit$len
    X <- fit$X
    estMat <- fit$estMat
    tsOrder <- fit$tsOrder
    #adjust scaled parameters back to original scale
    intercepts <- fit$intercepts
    covMeans <- matrix(0, nrow = p, ncol = tsOrder)
    covScaleFactors <- matrix(0, nrow = p, ncol = tsOrder)
    #store coefficient column vectors for each slice 1-d
    betasOrigScale <- array(0, c(p, p, tsOrder))
    intsOrigScale <- intercepts
    if (tsOrder >= 1)
    {
      for (j in 1:tsOrder)
      {
        means <- apply(X[,,len-j], 2, mean)
        scaleFactors <- apply(X[,,len-j], 2, sd)*sqrt((n-1)/n)
        coefsOrigScale <- t(estMat[,,(d-j+1)])/scaleFactors
        betasOrigScale[,,j] <- coefsOrigScale
        intsOrigScale <- intsOrigScale - apply(coefsOrigScale*means, 2, sum)
      }
    }

    i <- 0
    while (i <= tp)
    {
      Y <- matrix(rep(intsOrigScale, each = n), nrow = n, ncol = p)
      d1 <- dim(X)[3]
      if (tsOrder >= 1)
      {
        for (j in 1:tsOrder)
        {
          Y <- Y + X[,,(d1-j)]%*%betasOrigScale[,,j]
        }
      }
      X2 <- array(0, c(n, p, d1+1))
      X2[,,1:d1] <- X
      X2[,,(d1+1)] <- Y
      X <- X2
      rm(X2)
      i <- i+1
    }
    
    return(Y)
  }
