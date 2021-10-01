#' Creates the skeleton of a random sparse network
#' @param d lag of VAR process
#' @param p dimension of VAR process
#' @param n sample size
#' @param grp group structure of covariates
#' @return constructed d x p x p network
#' @export defn_net

defn_net =
  function(
    d,
    p,
    n,
    grp = NULL, 
    sparsity = NULL,
    grp_sparsity = 0.5
  ){
    if (is.null(grp))
    {
      grp <- 1:p
    }
    edge = array(0, c(d, p, p))
    weight = c(1, 1, 1)
    signum = c(1, -1)
    grpCt = length(unique(grp))
    
    if (is.null(sparsity)) {
      sparsity = max(min((n/(d*grpCt*p)), (0.05)), 0.01)
    }
    # cat(paste("sparsity =", round(sparsity, 4)))
    
    num_sel = ceiling(sparsity*(d*p*p))
    
    if (any(grp != 1:p)) {
      grp_sel = sample(1:grpCt, floor(grpCt * grp_sparsity))
      j_grp = which(grp %in% grp_sel)
      num_pergrp = floor(num_sel / length(j_grp) )
      
      for (j in j_grp){
        ind_sel = sample(((j-1)*p*d+1):(d*p*j), num_pergrp)
        edge[ind_sel] = sample(signum, num_pergrp, replace = TRUE) * 0.5
      }
    } else {
      ind_sel = sample(1:(p*p*d), num_sel)
      edge[ind_sel] = sample(signum, num_sel, replace = TRUE)* 0.5
    }
    
    bottom = cbind(diag(rep(1,p*(d-1))), matrix(0,p*(d-1),p))
    while(1) {
      A = NULL
      for (i in 1:d){
        A = cbind(A,Matrix(edge[i,,]))
      }
      At = rbind(A, bottom)
      if (max(abs(eigen(At)$values)) < 0.95) {
        # absolute summable, y_t is var(1) process
        break
      }
      edge = edge*0.95
    }
    return(edge)
  }


