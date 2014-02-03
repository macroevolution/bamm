#############################################################
#
#	branchMeanRateExponential(....)
#
#	Computes time-averaged rate of eponential process

branchMeanRateExponential <- function(t1, t2, p1, p2){
	res <- vector(mode = 'numeric', length = length(t1));
	res[which(p2 == 0)] <- p1[which(p2 == 0)];
	nonzero <- which(p2 < 0);
	res[nonzero] <- (p1[nonzero]/p2[nonzero])*(exp(p2[nonzero]*t2[nonzero]) - exp(p2[nonzero]*t1[nonzero])) / (t2[nonzero] - t1[nonzero]);
	nonzero <- which(p2 > 0);
	res[nonzero] <- (p1[nonzero]/p2[nonzero])*(2*p2[nonzero]*(t2[nonzero]-t1[nonzero]) + exp(-p2[nonzero]*t2[nonzero]) - exp(-p2[nonzero]*t1[nonzero])) / (t2[nonzero] - t1[nonzero]);
	return(res);
}