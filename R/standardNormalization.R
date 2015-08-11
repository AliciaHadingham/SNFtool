# Normalize each column of x to have mean 0 and standard deviation 1.


#' Standard Normalization
#' 
#' Normalize each column of the input data to have mean 0 and standard
#' deviation 1.
#' 
#' 
#' @param x The unnormalized data.
#' @return The data normalized.
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
#' Data1 = standardNormalization(Data1);
#' Data2 = standardNormalization(Data2);
#' 
standardNormalization = function(x) {
  x = as.matrix(x);
  mean = apply(x, 2, mean)
  sd = apply(x, 2, sd)
  sd[sd==0] = 1
  xNorm = t((t(x) - mean) / sd)
  return(xNorm)
}
