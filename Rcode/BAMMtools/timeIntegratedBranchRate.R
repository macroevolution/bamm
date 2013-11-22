#############################################################
#
#	timeIntegratedBranchRate(....)
#		computes the integral of rates on a branch segment with respect to time
#		Not the average.
#		Integrates the exponential function p1 * exp(p2 * t)
#		

timeIntegratedBranchRate <- function(t1, t2, p1, p2){
	
	if (p2 == 0){
		return(p1 * (t2 - t1));
	}else{
		(p1/p2)*(exp(p2*t2) - exp(p2*t1));
	}
}
timeIntegratedBranchRate <- Vectorize(timeIntegratedBranchRate);
