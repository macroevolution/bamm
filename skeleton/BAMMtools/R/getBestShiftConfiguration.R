
getBestShiftConfiguration <- function(x, threshold=0.01){
	

	if (class(x) == 'bammdata'){
		x <- credibleShiftSet(x, threshold);	
	}else if (class(x) == 'credibleshiftset'){

	}else{
		stop("Argument x must be of class bammdata or credibleshiftset\n");
	}
	class(x) <- 'bammdata';	
	subb <- subsetEventData(x, index=1);
	
	return(subb);
 
}


