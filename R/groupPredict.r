#' Group Predict
#' 
#' This function is used to predict the subtype of new patients.
#' 
#' 
#' @param train Training data. Has the same number of view and columns as test
#' data.
#' @param test Test data. Has the same number of view and columns as training
#' data.
#' @param groups The label for the training data.
#' @param K Number of neighbors.
#' @param alpha Hyperparameter used in constructing similarity network.
#' @param t Number of iterations.
#' @param method A indicator of which method to use to predict the label.
#' method = 0 means to use local and global consistency; method = 1 means to
#' use label propagation.
#' @return Returns the prediction of which group the test data belongs to.
#' @author Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
#' @examples
#' 
#' 
#' # Provide an example of predicting the new labels with label propagation
#' 
#' # Load views into list "dataL" and the cluster assignment into vector "label"
#' data(dataL)
#' data(label)
#' 
#' # Create the training and test data
#' n = floor(0.8*length(label)) # number of training cases
#' trainSample = sample.int(length(label), n)
#' train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
#' test = lapply(dataL, function(x) x[-trainSample, ]) # Test the rest of the data set
#' groups = label[trainSample]
#' 
#' # Set the other
#' K = 20
#' alpha = 0.5
#' t = 20
#' method = TRUE
#' 
#' # Apply the prediction function to the data
#' newLabel = groupPredict(train,test,groups,K,alpha,t,method)
#' 
#' # Compare the prediction accuracy
#' accuracy = sum(label[-trainSample] == newLabel[-c(1:n)])/(length(label) - n)
#' 
#' 
groupPredict <- function(train,test,groups,K=20,alpha=0.5,t=20,method=1){

###This function is used to predict the subtype of new patients.
#train and test have the same number of view and the same number of columns
# group is the label for the train data
# K, alpha, t are the prameters for SNF. 
#K is the number of neighbors
#alpha is the hyperparameter used in constructing similarity network
# t is the number of iterations
#method is a indicator of which method to use to predict the label. method = 0 means to use local and global consistency; method = 1 means to use label propagation.
Wi= vector("list", length=length(train));

for (i in 1:length(train)){
view= standardNormalization(rbind(train[[i]],test[[i]]));
Dist1 = dist2(view, view);
Wi[[i]] = affinityMatrix(Dist1, K, alpha);
}

W = SNF(Wi,K,t);
Y0=matrix(0,nrow(view), max(groups));
for (i in 1:length(groups)) Y0[i,groups[i]]=1;
Y=.csPrediction(W,Y0,method);
newgroups=rep(0,nrow(view));
for (i in 1:nrow(Y)) newgroups[i]=which(Y[i,]==max(Y[i,]));

return (newgroups);
}
