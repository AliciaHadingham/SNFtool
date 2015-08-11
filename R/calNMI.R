#' Mutual Information calculation
#' 
#' Calculate the mutual information between vectors x and y.
#' 
#' 
#' @param x a vector
#' @param y a vector
#' @return Returns the mutual information between vectors x and y.
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @references B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B
#' Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and
#' effective method to aggregate multiple data types on a genome wide scale.
#' Nature Methods. Online. Jan 26, 2014
#' @examples
#' 
#' 
#' # How to use SNF with multiple views
#' 
#' # Load views into list "dataL"
#' data(dataL)
#' data(label)
#' 
#' # Set the other parameters
#' K = 20 # number of neighbours
#' alpha = 0.5 # hyperparameter in affinityMatrix
#' T = 20 # number of iterations of SNF
#' 
#' # Normalize the features in each of the views if necessary
#' # dataL = lapply(dataL, standardNormalization)
#' 
#' # Calculate the distances for each view
#' distL = lapply(dataL, function(x) dist2(x, x))
#' 
#' # Construct the similarity graphs
#' affinityL = lapply(distL, function(x) affinityMatrix(x, K, alpha))
#' 
#' # Example of how to use SNF to perform subtyping
#' # Construct the fused network
#' W = SNF(affinityL, K, T)
#' # Perform clustering on the fused network.
#' clustering = spectralClustering(W,3);
#' # Use NMI to measure the goodness of the obtained labels.
#' NMI = calNMI(clustering,label);
#' 
#' 
calNMI <- function(x, y) {
##This function calculate the NMI between two clusters.	
	
	x = as.vector(x);
	y = as.vector(y);  
    return(max(0, .mutualInformation(x, y)/sqrt(.entropy(x) * .entropy(y)), na.rm=TRUE))
}


