#' Quantize values
#' 
#' Replaces values by their quantile number.
#' 
#' 
#' @param x numeric vector to quantize
#' @param quantiles number of quantiles to use
#' @return a numeric vector of the same length as x, giving the quantiles of values in x
#' @examples
#' 
#' quantize(runif(100))
#' 
#' # ranking
#' quantize(runif(100), 100)
#' 
#' # quantization can deal with high kurtosis data
#' quantize(c(rnorm(100,-1000), rnorm(100,1000)))
#' 
quantize <- function(x, quantiles=20){
  table <- quantile(x,seq(0,1,length.out=quantiles+1),names=FALSE, na.rm=TRUE)
  findInterval(x,sort(table),all.inside=TRUE)
}



#' Create a similarity matrix, robustly
#' 
#' This function uses quantization to calculate useful similarity matricies, even with poorly behaved data.
#' 
#' 
#' @param x a numeric data frame, array or matrix. Rows are samples, columns are attributes
#' @return a similarity matrix
#' @examples
#' 
#' robustSimilarity(mtcars)
#' robustSimilarity(as.matrix(mtcars))
#' 
robustSimilarity <- function(x){
  quantiles <- min(nrow(x), 20)
  distances <- as.matrix(dist(apply(x,2,quantize,quantiles=quantiles)))
  similarities <- exp(-((distances/mean(distances))^2))
  return(similarities)
}
