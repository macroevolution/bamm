getEventCorrelationMatrix <-
function(ephy){

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
