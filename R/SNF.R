SNF <- function(Wall,K=20,t=20,layer_bias = rep.int(1, length(Wall)),parallel=FALSE, auto_stop=FALSE) {
	
	###This function is the main function of our software. The inputs are as follows:
    # Wall : List of affinity matrices
    # K : number of neighbors
    # t : number of iterations for fusion
    # layer_bias : numerical vector of length length(Wall) containing the weights for each layer
    # parallel : should the algorithm attempt to run in parallel? A parallel backend needs to be set up first?
    # auto_stop : should the algorithm stop early if it believes it has converged?
    
    ###The output is a unified similarity graph. It contains both complementary information and common structures from all individual network. 
    ###You can do various applications on this graph, such as clustering(subtyping), classification, prediction.
    
    LW = length(Wall)
    normalize <- function(X) X / rowSums(X)
    # makes elements other than largest K zero
    
    
    newW <- vector("list", LW)
    nextW <- vector("list", LW)
    ###First, normalize different networks to avoid scale problems.
    for( i in 1: LW){
      Wall[[i]] = normalize(Wall[[i]]);
      Wall[[i]] = (Wall[[i]]+t(Wall[[i]]))/2;
    }
    
    ### Calculate the local transition matrix.
    for( i in 1: LW){
      newW[[i]] = (.dominateset(Wall[[i]],K))
    }
    
    # Set up convergence detection
    devs <- NULL
    converged <- FALSE
    
    # perform the diffusion for t iterations
    for (i in 1:t) {
      nextW <- plyr::llply(1:LW, function(j){
        sumWJ = matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
        adjusted_bias <- layer_bias/sum(layer_bias[-j])
        for (k in setdiff(1:LW, j)) {
          sumWJ = sumWJ + (Wall[[k]]*adjusted_bias[k])
        }
        return(newW[[j]] %*% (sumWJ/(LW - 1)) %*% t(newW[[j]]))
      }, .parallel=parallel)
      ###Normalize each new obtained networks.
      for(j in 1 : LW){
      	      
                Wall[[j]] = nextW[[j]] + diag(nrow(Wall[[j]]));
                Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2;
                }
      
      
      # Check for convergence
      wallarray <- array(unlist(Wall),dim=c(dim(Wall[[1]]),length(Wall)))
      means <- base::rowMeans(wallarray, dims=2)
      dev <- mean(abs(wallarray-array(means,dim=dim(wallarray))) / mean(means))
      devs <- c(devs,dev)
      d_devs <- (devs-lag(devs))/lag(devs)
      d2_devs <- (d_devs-lag(d_devs))/lag(d_devs)
      if(!is.na(d2_devs[length(d2_devs)]) & 
           abs(d_devs[length(d_devs)])<0.01 & 
           abs(d2_devs[length(d2_devs)])<0.01){
        converged <- TRUE
      }
      if(auto_stop & converged){
        break
      }
   }
    
    # construct the combined affinity matrix by summing diffused matrices
    W = matrix(0,nrow(Wall[[1]]), ncol(Wall[[1]]))
    for(i in 1:LW){
      W = W + Wall[[i]]
    }
    W = W/LW;
    W = normalize(W);
    # ensure affinity matrix is symmetrical
    W = (W + t(W)+diag(nrow(W))) / 2;
    
    return(W)  
  }
