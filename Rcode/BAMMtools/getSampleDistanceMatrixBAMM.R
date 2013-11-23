#############################################################
#
#	getSampleDistanceMatrixBAMM <- function(...)
#   For each of k samples from the posterior with branch-specific rates
#		computes pairwise rate differences between sample
#		Thus, for any two rate configurations i and j
#			the i,j and j,i elements of the matrix will be 
#			the squared rate differences
#	
#		phylist is a list of phylogenetic trees with branch-specific rates for 
#			branch lengths (rather than time)

getSampleDistanceMatrixBAMM <- function(phylist, modeltree){
	
	
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

