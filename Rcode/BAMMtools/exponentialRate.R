#############################################################
#
#	exponentialRate(....)
#
#	Gets point estimate of evolutionary rate
#	Vectorized

exponentialRate <- function(t1, p1, p2){
	(p1 * exp(p2 * t1));
}
exponentialRate <- Vectorize(exponentialRate);
