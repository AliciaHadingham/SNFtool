#' Estimate Number Of Clusters Given Graph
#' 
#' This function estimates the number of clusters given the two huristics given
#' in the supplementary materials of our nature method paper W is the
#' similarity graph NUMC is a vector which contains the possible choices of
#' number of clusters.
#' 
#' 
#' @param W List of matrices. Each element of the list is a square, symmetric
#' matrix that shows affinities of the data points from a certain view.
#' @param NUMC A vector which contains the possible choices of number of
#' clusters.
#' @return K1 is the estimated best number of clusters according to eigen-gaps
#' K12 is the estimated SECOND best number of clusters according to eigen-gaps
#' K2 is the estimated number of clusters according to rotation cost K22 is the
#' estimated SECOND number of clusters according to rotation cost
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @references B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B
#' Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and
#' effective method to aggregate multiple data types on a genome wide scale.
#' Nature Methods. Online. Jan 26, 2014
#' 
#' Concise description can be found here:
#' http://compbio.cs.toronto.edu/SNF/SNF/Software.html
#' @examples
#' 
#' 
#' ## First, set all the parameters:
#' K = 20;  	# number of neighbors, usually (10~30)
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
#' ## These similarity graphs have complementary information about clusters.
#' displayClusters(W1,truelabel);
#' displayClusters(W2,truelabel);
#' 
#' ## next, we fuse all the graphs
#' ## then the overall matrix can be computed by similarity network fusion(SNF):
#' W = SNF(list(W1,W2), K, T)
#' 
#' ## With this unified graph W of size n x n, 
#' ## you can do either spectral clustering or Kernel NMF. 
#' ## If you need help with further clustering, please let us know. 
#' 
#' ## You can display clusters in the data by the following function
#' ## where C is the number of clusters.
#' C = 2 								# number of clusters
#' group = spectralClustering(W,C); 	# the final subtypes information
#' displayClusters(W, group)
#' 
#' ## You can get cluster labels for each data point by spectral clustering
#' labels = spectralClustering(W, C)
#' 
#' plot(Data1, col=labels, main='Data type 1')
#' plot(Data2, col=labels, main='Data type 2')
#' 
#' ## Here we provide two ways to estimate the number of clusters. Note that,
#' ## these two methods cannot guarantee the accuracy of esstimated number of
#' ## clusters, but just to offer two insights about the datasets.
#' 
#' estimationResult = estimateNumberOfClustersGivenGraph(W, 2:5);
#' 
estimateNumberOfClustersGivenGraph <- function(W, NUMC=2:5) {
  
  #   This function estimates the number of clusters given the two huristics
  #   given in the supplementary materials of our nature method paper
  #   W is the similarity graph
  #   NUMC is a vector which contains the possible choices of number of
  #   clusters.
  #   
  #   
  #   K1 is the estimated best number of clusters according to eigen-gaps
  #   K12 is the estimated SECOND best number of clusters according to eigen-gaps
  #   
  #   K2 is the estimated number of clusters according to rotation cost
  #   K22 is the estimated SECOND number of clusters according to rotation cost
  #   
  #   an example would be [K1, K2, K12,K22] = Estimate_Number_of_Clusters_given_graph(W,
  #                                                                                   [2:5]);
  #   
  #   Note that this function can only give an estimate of the number of
  #   clusters. How to determine the "OPTIMAL" number of clusters, is still an
  #   open question so far. 
  
  if (min(NUMC) == 1) {
    warning('Note that we always assume there are more than one cluster.');
    NUMC = NUMC[NUMC > 1]  
  }
  
  W = (W + t(W))/2
  diag(W) = 0
  
  if (length(NUMC) > 0) {
    degs = rowSums(W)

    
    # compute unnormalized Laplacian
    
    degs[degs == 0] = .Machine$double.eps    
    D = diag(degs)    
    L = D - W
    Di = diag(1 / sqrt(degs))
    L = Di %*% L %*% Di
    
    # compute the eigenvectors corresponding to the k smallest
    # eigs$valuess
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return=T)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = eigengap * (1 - eigs$values[1:length(eigs$values) - 1] ) / (1 - eigs$values[2:length(eigs$values)])

    quality = list()
    for (c_index in 1:length(NUMC)) {
      ck = NUMC[c_index]
      UU = eigs$vectors[, 1:ck]
      EigenvectorsDiscrete <- .discretisation(UU)[[1]]
      EigenVectors = EigenvectorsDiscrete^2
      
      # MATLAB: sort(EigenVectors,2, 'descend');
      temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), function(i) EigenVectors[, i])), ]
      temp1 <- t(apply(temp1, 1, sort, TRUE))  
      
      quality[[c_index]] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) * 
        sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*% temp1[, 1:max(2, ck-1)] ))
    }
    
    t1 <- sort(eigengap[NUMC], decreasing=TRUE, index.return=T)$ix
    K1 = NUMC[t1[1]]
    K12 = NUMC[t1[2]]
    t2 <- sort(unlist(quality), index.return=TRUE)$ix
    K2 <- NUMC[t2[1]]
    K22 <- NUMC[t2[2]]    
  }
  
  return (list(K1, K12, K2, K22))
}
