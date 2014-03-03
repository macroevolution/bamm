
getBestShiftConfiguration <- function(x, threshold){
	
	if (class(threshold) != 'branchprior' & class(threshold) != 'numeric'){
		stop("arg threshold is not valid\n")
	}
	
	if (class(x) == 'bammdata'){
		x <- credibleShiftSet(x, set.limit = 0.95, threshold = threshold);	
	}else if (class(x) == 'credibleshiftset'){

	}else{
		stop("Argument x must be of class bammdata or credibleshiftset\n");
	}
	class(x) <- 'bammdata';	
	subb <- subsetEventData(x, index=1);
	
	return(subb);
 
}


