################################################################################
#' @description convert gene expression matrix to rank matrix denoted by R with element r_is, which represents the rank of gene g_i in sample s
#' @param X Expression matrix data of samples.
#' @return rank matrix
################################################################################
rank.matrix <- function(x){
  rankmatrix = c()
  for (i in 1: ncol(x)){
    temp = rank(x[,i])
    rankmatrix = cbind(rankmatrix, temp)
  }
  colnames(rankmatrix) = colnames(x)
  row.names(rankmatrix) = row.names(x)
  return(rankmatrix)
}	



#####################################################################
#' @description Remove the Null element in a list
#' @param X A list with Null elements
#' @return A list without Null elements
######################################################################
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}



########################################################################################
#' @description caculate the Delta Rank matrix for all samples
#' @param X=[X1,X2,X3], x1 is the rank matrix for normal samles, x2 is the rank vector for mean value of X1£¬x3 is the rank matrix for cancer samles
#' @param net the backgound network (two column genes assocaited in the background network)
#' @param n_normal the number of normal samples,i.e. ncol(X1)
#' @param n_cancer the number of cancer samples,i.e. ncol(X3), ncol(X) = n_noraml + 1 + n_cancer
#' @return net.edge, edges of background network;
#' @return net.data, deltarank vector for mean of normal samples
#' @return deltarank_normal,  deltarank matrix for normal samples
#' @return deltarank_cancer,  deltarank matrix for cancer samples
#########################################################################################
delta.rank <- function(net, x, n_normal, n_cancer){
  deltarank = c()
  deltarank <- sapply(1: nrow(net), function(i){
    print(i)
    r1= which(row.names(x) == net[i,1])
    r2= which(row.names(x) == net[i,2])
    if((length(r1)!= 0) & (length(r2)!= 0)){
      r <- x[r1,]-x[r2,]
      deltarank <- rbind(deltarank, c(as.matrix(net[i,]), r))
    }
  })
  deltarank <- rmNullObs(deltarank)
  deltarank<- do.call(rbind, deltarank)
  net <- matrix(unlist(deltarank[, c(1, 2, n_normal+3)]), ncol=3) 
  net.edge = matrix(unlist(net[,c(1,2)]), ncol=2)
  net.data = matrix(unlist(net[,3]), ncol=1)
  deltarank_normal = matrix(unlist(deltarank[, 3:(n_normal+2)]), ncol = n_normal)
  deltarank_cancer = matrix(unlist(deltarank[, (n_normal+4):ncol(deltarank)]), ncol=n_cancer)
  net.data = apply(net.data, 2, as.numeric)
  deltarank_normal = apply(deltarank_normal, 2, as.numeric)
  colnames(deltarank_normal)  = colnames(x)[1:n_normal]
  deltarank_cancer = apply(deltarank_cancer, 2, as.numeric)
  colnames(deltarank_cancer)  = colnames(x)[(n_normal+2): ncol(x)]
  return(list(net.edge, net.data, deltarank_normal, deltarank_cancer))
}

########################################################################################
#' @description caculate the edge-perturbation matrix (EPm)
#' @param deltarank deltarank matrix; net.data, deltarank vector for mean of normal samples
#' @return EPm, specific delta Rank matrix
#########################################################################################
EPm <- function(net.data, deltarank){
  EPm <- c()
  EPm <- sapply( 1: ncol(deltarank), function(i){
    delta <- deltarank[,i] - stablenet.data
    EPm <- cbind(EPm, delta)
  })
  colnames(EPm)  <- colnames(deltarank)
  return(EPm)
}






