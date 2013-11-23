getSampleCoMat <- function(phylist, modeltree){
	
	comat <- matrix(NA, nrow=length(phylist[[1]]$edge.length), ncol=length(phylist));
	
	# first, scale all edge lengths by time
	for (i in 1:length(phylist)){
		comat[,i] <- phylist[[i]]$edge.length * modeltree$edge.length;
	}	
 
	return(comat);
	
}
