#' Affinity matrix calculation
#' 
#' Computes affinity matrix from a generic distance matrix
#' 
#' 
#' @param Diff Distance matrix
#' @param K Number of nearest neighbors
#' @param sigma Variance for local model
#' @return Returns an affinity matrix that represents the neighborhood graph of
#' the data points.
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @references B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B
#' Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and
#' effective method to aggregate multiple data types on a genome wide scale.
#' Nature Methods. Online. Jan 26, 2014
#' @examples
#' 
#' 
#' ## First, set all the parameters:
#' K = 20; ##number of neighbors, must be greater than 1. usually (10~30)
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
#' ## Calculate distance matrices(here we calculate Euclidean Distance, 
#' ## you can use other distance, e.g. correlation)
#' Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
#' Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))
#' 
#' ## Next, construct similarity graphs
#' W1 = affinityMatrix(Dist1, K, alpha)
#' W2 = affinityMatrix(Dist2, K, alpha)
#' 
#' 
affinityMatrix <- function(Diff,K=20,sigma=0.5) {
###This function constructs similarity networks.
  N = nrow(Diff)
  
  Diff = (Diff + t(Diff)) / 2
  diag(Diff) = 0;
  sortedColumns = as.matrix(t(apply(Diff,2,sort)))
  finiteMean <- function(x) { mean(x[is.finite(x)]) }
  means = apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;
  
  avg <- function(x,y) ((x+y)/2)
  Sig = outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
  Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
  densities = dnorm(Diff,0,sigma*Sig,log = FALSE)
  
  W = (densities + t(densities)) / 2
  return(W)
}
