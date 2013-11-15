getBySpeciesRateMatrix <-
function(ephy, start.time, nbreaks, ndr=TRUE, node){
	
	
	spset <- ephy$tip.label;
	
	mm <- matrix(NA, nrow=length(spset), ncol=nbreaks);
	for (i in 1:nrow(mm)){
		#cat(spset[i], '\n')
		mm[i,] <- getSpeciesRateThroughTime(ephy, start.time, nbreaks, ndr, species=spset[i]);
	}
	
	rownames(mm) <- spset;
	return(mm);
}
