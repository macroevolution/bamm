#############################################################
#
#	getEventCorrelationMatrix(....)
#
# 	Each entry of this matrix represents the expected 
#	probability that a pair[i, j] of tips will have the same 
#	rate parameters due to BAMM model.
#
# 	Should modify this to allow exponential, spherical,
#		and other possible correlation structures.
#	Need to make a corStruct class that works with this
#		for GLS analyses

getEventCorrelationMatrix <- function(ephy){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	TOL <- 0.0001;
	corMat <- matrix(0, nrow=length(ephy$tip.label), ncol=length(ephy$tip.label));
	for (i in 1:length(ephy$tipStates)){
		dd <- dist(ephy$tipStates[[i]]);
		cmat <- as.matrix(dd);	
		corMat <- corMat + (cmat < TOL);
		
	}
	rownames(corMat) <- ephy$tip.label;
	colnames(corMat) <- ephy$tip.label;
	return(corMat/ length(ephy$numberEvents));
}
