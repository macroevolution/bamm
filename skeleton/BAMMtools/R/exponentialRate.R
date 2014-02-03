#############################################################
#
#	exponentialRate(....)
#
#	Gets point estimate of evolutionary rate
#	Vectorized <- No. this turns out to be much much slower

exponentialRate <- function(t1, p1, p2) {
	zero <- which(p2 == 0);
	ret <- numeric(length(t1));
	ret[zero] <- p1[zero];
	nonzero <- which(p2 != 0);
	p1 <- p1[nonzero];
	p2 <- p2[nonzero];
	t1 <- t1[nonzero];
	ret[nonzero] <- ( p1 * ((p2/abs(p2)) * (1 - exp(-abs(p2)*t1)) + 1) );
	return(ret);
}
#exponentialRateV <- Vectorize(exponentialRate);
