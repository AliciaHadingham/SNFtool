## Arguments:
## data: a list, where each item in the list is a matrix of values for each data type
## W: the target network for which the NMI is calculated against for each feature
##
## Details:
## NMI is calculated based on the clustering assignments using spectral clustering
## The number of clusters is set based on the estimateNumberOfClustersGivenGraph on the target matrix 
## using default parameters.
##
## Outputs:
## A list that contains the NMI score for each feature and their ranks from highest to lowest
## output[[1]] is the NMI score
## output[[1]][[1]] is the NMI score of first data type
## output[[1]][[1]][1] is the NMI score of the first feature of the first data type
## similarly for output[[2]]... except it is the rank instead of the score



#' Rank Features by NMI
#' 
#' Ranks each features by NMI based on their clustering assingments
#' 
#' 
#' @param data List containing all the data types.
#' @param W Target Matrix for which the NMI is calculated against.
#' @return List containing the NMI and rank based on NMI for each feature.
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @examples
#' 
#' 
#' ## First, set all the parameters:
#' K = 20;		# number of neighbors, usually (10~30)
#' alpha = 0.5;  	# hyperparameter, usually (0.3~0.8)
#' T = 20; 	# Number of Iterations, usually (10~20)
#' 
#' ## Data1 is of size n x d_1, 
#' ## where n is the number of patients, d_1 is the number of genes, 
#' ## Data2 is of size n x d_2, 
#' ## where n is the number of patients, d_2 is the number of methylation
#' data(Data1)
#' data(Data2)
#' 
#' ## Here, the simulation data (SNFdata) has two data types. They are complementary to each other. 
#' ## And two data types have the same number of points. 
#' ## The first half data belongs to the first cluster; the rest belongs to the second cluster.
#' truelabel = c(matrix(1,100,1),matrix(2,100,1)); ## the ground truth of the simulated data
#' 
#' ## Calculate distance matrices
#' ## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
#' 
#' ## If the data are all continuous values, we recommend the users to perform 
#' ## standard normalization before using SNF, 
#' ## though it is optional depending on the data the users want to use.  
#' # Data1 = standardNormalization(Data1);
#' # Data2 = standardNormalization(Data2);
#' 
#' 
#' 
#' ## Calculate the pair-wise distance; 
#' ## If the data is continuous, we recommend to use the function "dist2" as follows 
#' Dist1 = dist2(as.matrix(Data1),as.matrix(Data1));
#' Dist2 = dist2(as.matrix(Data2),as.matrix(Data2));
#' 
#' ## next, construct similarity graphs
#' W1 = affinityMatrix(Dist1, K, alpha)
#' W2 = affinityMatrix(Dist2, K, alpha)
#' 
#' ## next, we fuse all the graphs
#' ## then the overall matrix can be computed by similarity network fusion(SNF):
#' W = SNF(list(W1,W2), K, T)
#' 
#' NMI_scores <- rankFeaturesByNMI(list(Data1, Data2), W)
#' 
#' 
rankFeaturesByNMI <- function(data, W) 
{  
  stopifnot(class(data) == "list")
  
  NUM_OF_DATA_TYES <- length(data)
  NMI_scores <- vector(mode="list", length=NUM_OF_DATA_TYES)
  NMI_ranks <- vector(mode="list", length=NUM_OF_DATA_TYES)
  num_of_clusters_fused <- estimateNumberOfClustersGivenGraph(W)[[1]]
  clustering_fused <- spectralClustering(W, num_of_clusters_fused)
  
  for (data_type_ind in 1:NUM_OF_DATA_TYES)
  {
    NUM_OF_FEATURES <- dim(data[[data_type_ind]])[2] 
    NMI_scores[[data_type_ind]] <- vector(mode="numeric", length=NUM_OF_FEATURES)    
    
    for (feature_ind in 1:NUM_OF_FEATURES)
    {
      affinity_matrix <- affinityMatrix(
        dist2(as.matrix(data[[data_type_ind]][, feature_ind]), as.matrix(data[[data_type_ind]][, feature_ind])))      
      clustering_single_feature <- spectralClustering(affinity_matrix, num_of_clusters_fused)
      NMI_scores[[data_type_ind]][feature_ind] <- calNMI(clustering_fused, clustering_single_feature)      
    }    
    NMI_ranks[[data_type_ind]] <- rank(-NMI_scores[[data_type_ind]], ties.method="first")
  }
  
  return(list(NMI_scores, NMI_ranks))  
}
