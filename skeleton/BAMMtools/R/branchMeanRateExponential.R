#############################################################
#
#	branchMeanRateExponential(....)
#
#	Computes time-averaged rate of eponential process
#	Vectorized


branchMeanRateExponential <- function(t1, t2, p1, p2){
	if (p2 == 0){
		return(p1);
	}else{
		(p1/p2)*(exp(p2*t2) - exp(p2*t1)) / (t2 - t1);
	}
}
branchMeanRateExponential <- Vectorize(branchMeanRateExponential);
