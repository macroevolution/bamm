safeCompare.old = function(vec,val,FUN,tol=1e-4) {
	if (FUN == ">=") {
		ret = rep(FALSE,length(vec));
		ret[(vec-val) >= 0] = TRUE;
		ret[abs(vec-val) <= tol] = TRUE;
		return(ret);
	}
	if (FUN == "<=") {
		ret = rep(FALSE,length(vec));
		ret[(vec-val) <= 0] = TRUE;
		ret[(vec-val) <= tol] = TRUE;
		return(ret);
	}
}


safeCompare <- function (vec, val, comp, tol=1e-4) {
	if (comp == ">=") {
		return((val-vec) <= tol);
	} else if (comp == "<=") {
		return((vec-val) <= tol);
	} else {
		stop("Invalid 'comp'.");
	}
}