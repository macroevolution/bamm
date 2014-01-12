#############################################################
#
#	branchMeanRateExponential(....)
#
#	Computes time-averaged rate of eponential process

branchMeanRateExponential <- function(t1, t2, p1, p2){
	res <- vector(mode = 'numeric', length = length(t1));
	res[which(p2 == 0)] <- p1[which(p2 == 0)];
	nonzero <- which(p2 != 0);
	p1 <- p1[nonzero];
	p2 <- p2[nonzero];
	t1 <- t1[nonzero];
	t2 <- t2[nonzero];
	res[nonzero] <- (p1/p2)*(exp(p2*t2) - exp(p2*t1)) / (t2 - t1);
	return(res);
}