
getSampleCorrelationMatrix <- function(ephy, verbose=F) {

	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}
	
	TOL <- 0.0001;
	corMat <- matrix(0, nrow=length(ephy$eventData), ncol=length(ephy$eventData));

	### Step 1: make a list of within-sample distance matrices
	
	xx <- list();
	for (i in 1:length(ephy$eventData)){
		
		tmp <- as.matrix(dist(ephy$tipStates[[i]])); 
		xx[[i]] <- 1*(tmp < TOL);
	}

	# Step 2: Compute matrix correlation between samples i and j
	# This really needs to be in C++ : very slow
	for (i in 1:(length(ephy$eventData)-1)){
		if (verbose == T){
			cat("Processing event ", i, '\n');		
		}

		for (j in (i + 1):length(ephy$eventData)){
			corMat[i,j] <- matrixCorr(xx[[i]], xx[[j]]);
			corMat[j, i] <- corMat[i,j];
			
		}
	}
	diag(corMat) <- rep(1, length(ephy$eventData));

	return(corMat);
}



