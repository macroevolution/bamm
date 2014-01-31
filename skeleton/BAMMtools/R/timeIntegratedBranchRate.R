#############################################################
#
#	timeIntegratedBranchRate(....)
#		computes the integral of rates on a branch segment with respect to time
#		Not the average.
#		Integrates the exponential function p1 * exp(p2 * t)
#		

timeIntegratedBranchRate <- function(t1, t2, p1, p2){
	res <- vector(mode = 'numeric', length = length(t1));
	zero <- which(p2 == 0);
	res[zero] <- p1[zero] * (t2[zero] - t1[zero]);
	nonzero <- which(p2 < 0);
	res[nonzero] <- (p1[nonzero]/p2[nonzero])*(exp(p2[nonzero]*t2[nonzero]) - exp(p2[nonzero]*t1[nonzero]));
	nonzero <- which(p2 > 0);
	res[nonzero] <- (p1[nonzero]/p2[nonzero])*(2*p2[nonzero]*(t2[nonzero]-t1[nonzero]) + exp(-p2[nonzero]*t2[nonzero]) - exp(-p2[nonzero]*t1[nonzero]));
	return(res);
}
