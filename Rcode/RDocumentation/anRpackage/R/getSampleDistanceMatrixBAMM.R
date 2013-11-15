getSampleDistanceMatrixBAMM <-
function(phylist, modeltree){
	
	
	# first, scale all edge lengths by time
	for (i in 1:length(phylist)){
		phylist[[i]]$edge.length <- phylist[[i]]$edge.length * modeltree$edge.length;
	}
	
	dmat <- matrix(0, nrow=length(phylist), ncol=length(phylist));
		
	nsamples <- length(phylist);	
	for (i in 1:(nsamples - 1)){
		for (j in i:nsamples){
			dmat[i, j] <- sum((phylist[[i]]$edge.length - phylist[[j]]$edge.length)^2);
		}
	}	
	dmat <- dmat + t(dmat);
	
	return(dmat);
	
}
