#' Spectral Clustering
#' 
#' Perform the famous spectral clustering algorithms. There are three variants.
#' The default one is the third type.
#' 
#' 
#' @param affinity Similarity matrix
#' @param K Number of clusters
#' @param type The variants of spectral clustering to use.
#' @return A vector consisting of cluster labels of each sample.
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @examples
#' 
#' 
#' ## First, set all the parameters:
#' K = 20;##number of neighbors, usually (10~30)
#' alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
#' T = 20; ###Number of Iterations, usually (10~50)
#' 
#' ## Data1 is of size n x d_1, 
#' ## where n is the number of patients, d_1 is the number of genes, 
#' ## Data2 is of size n x d_2, 
#' ## where n is the number of patients, d_2 is the number of methylation
#' data(Data1)
#' data(Data2)
#' 
#' ## Calculate distance matrices (here we calculate Euclidean Distance, 
#' ## you can use other distance, e.g. correlation)
#' Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
#' Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))
#' 
#' ## Next, construct similarity graphs
#' W1 = affinityMatrix(Dist1, K, alpha)
#' W2 = affinityMatrix(Dist2, K, alpha)
#' 
#' # Next, we fuse all the graphs
#' # then the overall matrix can be computed by
#' W = SNF(list(W1,W2), K, T)
#' 
#' ## With this unified graph W of size n x n, 
#' ## you can do either spectral clustering or Kernel NMF. 
#' ## If you need help with further clustering, please let us know. 
#' 
#' ## You can display clusters in the data by the following function
#' ## where C is the number of clusters.
#' C = 2
#' 
#' ## You can get cluster labels for each data point by spectral clustering
#' labels = spectralClustering(W, C)
#' 
spectralClustering <- function(affinity, K, type=3) {
	
  ###This function implements the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
  ###THe inputs are as follows:
  
      #affinity: the similarity matrix;
      #K: the number of clusters
      # type: indicators of variants of spectral clustering 
	
  d = rowSums(affinity)
  d[d == 0] = .Machine$double.eps
  D = diag(d)
  L = D - affinity
  if (type == 1) {
    NL = L
  } else if (type == 2) {
    Di = diag(1 / d)
    NL = Di %*% L
  } else if(type == 3) {
    Di = diag(1 / sqrt(d))
    NL = Di %*% L %*% Di
  }
  eig = eigen(NL)
  res = sort(abs(eig$values),index.return = TRUE)
  U = eig$vectors[,res$ix[1:K]]
  normalize <- function(x) x / sqrt(sum(x^2))
  if (type == 3) {
    U = t(apply(U,1,normalize))
  }
  eigDiscrete = .discretisation(U)
  eigDiscrete = eigDiscrete$discrete
  labels = apply(eigDiscrete,1,which.max)
  
  
 
  return(labels)
}
