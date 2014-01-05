safeCompare <-
function (vec, val, comp, tol=1e-4) {
	if (comp == ">=") {
		return((val-vec) <= tol);
	} else if (comp == "<=") {
		return((vec-val) <= tol);
	} else {
		stop("Invalid 'comp'.");
	}
}
