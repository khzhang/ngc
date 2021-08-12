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
    sparsity = NULL
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
    cat(paste("sparsity =", round(sparsity, 4)))
    
    for (ii in 1:d){
      for (i in 1:p){
        for (j in unique(grp)){
            edge[ii, i, grp == j] = ((runif(1, 0, 1)< sparsity))*sample(weight, 1)*sample(signum, 1)
        }
      }
    }
    return(edge)
  }
