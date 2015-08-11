#' Pairwise Chi-squared distances
#' 
#' Computes the Chi-squared distances between all pairs of data point given
#' 
#' 
#' @param A A data matrix where each row is a different data point
#' @param B A data matrix where each row is a different data point. If this
#' matrix is the same as X, pairwise distances for all data points are
#' computed.
#' @return Returns an N x M matrix where N is the number of rows in X and M is
#' the number of rows in M. element (n,m) is the squared Chi-squared distance
#' between nth data point in X and mth data point in C
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @examples
#' 
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
#' Dist1 = chiDist2(as.matrix(Data1), as.matrix(Data1))
#' Dist2 = chiDist2(as.matrix(Data2), as.matrix(Data2))
#' 
#' 
chiDist2 <- function(A,B){
	###This function implements the Chi-Square distance between A and B
  n = nrow(A)
  m = nrow(B)
   d = ncol(A)
   stopifnot(d == ncol(B))
   res = matrix(nrow = n, ncol = m)
   sqA = A^2
   sqB = B^2
   twoAB = 2 * (A %*% t(B))
   for (a_num in 1:n) {
   	for (b_num in 1:m) {
   		res[a_num, b_num] = sum((sqA[a_num, ] + sqB[b_num, ] - twoAB[a_num, b_num]) / (A[a_num, ] + B[b_num, ]))
   		}
   		}
   		res = res / 2
   		return(res)  
}
