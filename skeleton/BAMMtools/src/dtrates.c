#include <R.h>
#include <Rinternals.h>

//forward function declarations
SEXP getListElement(SEXP list, char *str);
SEXP dtrates(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample, SEXP type);
double getDblMatrixELT(SEXP matrix, int row, int col);
double getMeanRateExponential(double t1, double t2, double p1, double p2);
double getTimeIntegratedBranchRate(double t1, double t2, double p1, double p2);

SEXP getListElement(SEXP list, char *str) {
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < LENGTH(list); i++) {
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

double getDblMatrixELT(SEXP matrix, int row, int col) {
	int nrow = INTEGER(getAttrib(matrix, R_DimSymbol))[0];
	return REAL(matrix)[row + nrow*col];
}

double getMeanRateExponential(double t1, double t2, double p1, double p2) {	
	if (p2 == 0.) {
		return p1;
	}
	else if (p2 < 0) {
		return (p1/p2)*(exp(t2*p2) - exp(t1*p2))/(t2 - t1);
	}
	else {
		return (p1/p2)*(2.*p2*(t2-t1) + exp(-t2*p2) - exp(-t1*p2))/(t2 - t1);
	}
}

double getTimeIntegratedBranchRate(double t1, double t2, double p1, double p2) {
	if (p2 == 0.) {
		return (t2 - t1) * p1;
	}
	else if (p2 < 0) {
		return (p1/p2)*(exp(t2*p2) - exp(t1*p2));
	}
	else {
		return (p1/p2)*(2.*p2*(t2-t1) + exp(-t2*p2) - exp(-t1*p2));
	}
}

/***
ephy is the bammdata object, a list holding all the relevant data.
segmat is a matrix where the rows represent branches of a phylogeny
broken up into many small segments.  each row indexes one such segment.
columns 2 and 3 give the starting and ending times of that segment and
column 1 is the node of the phylogeny to which that segement belongs.
tol is a precision parameter used for comparing starting and ending 
times of approximating segments and starting and ending times of branches
or branch segments on the phylogeny.
***/


SEXP dtrates(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample, SEXP type) {
	double eps = REAL(tol)[0];
	
	int k, j, l, nprotect = 0;
	int nsamples = LENGTH(sample);	
	int nsegs = INTEGER(getAttrib(segmat, R_DimSymbol))[0];
	
	SEXP rates, erates;
	PROTECT(rates = allocVector(REALSXP, nsegs)); nprotect++;
	for(k = 0; k < nsegs; k++) {
		REAL(rates)[k] = 0.;
	}
	if (INTEGER(type)[0] == 0) {
		PROTECT(erates = allocVector(REALSXP, nsegs)); nprotect++;
		for(k = 0; k < nsegs; k++) {
			REAL(erates)[k] = 0.;
		}	
	}
	
	int nrow, node, nnode, event, nxtevent, isGoodStart, isGoodEnd, place_holder;
	double begin, end, Start, lam1, lam2, mu1, mu2, relStart, relEnd, rightshift, leftshift, erightshift, eleftshift, ret;		
	
	for (k = INTEGER(sample)[0] - 1; k < INTEGER(sample)[nsamples - 1]; k++) {
		SEXP eventSegs, eventData;
		
		eventSegs = PROTECT(VECTOR_ELT(getListElement(ephy, "eventBranchSegs"), k)); nprotect++;
		eventData = PROTECT(VECTOR_ELT(getListElement(ephy, "eventData"), k)); nprotect++;
				
		nrow = INTEGER(getAttrib(eventSegs, R_DimSymbol))[0];
		place_holder = 0;
		//move down the rows of eventSegs
		for(j = 0; j < nrow; j++) {
			//eventSegs is 4 column matrix, node is in first column stored as double
			node = (int) REAL(eventSegs)[j + nrow * 0];
			
			//event index is in fourth column stored as double
			event = (int) REAL(eventSegs)[j  + nrow * 3];
			
			//begin and end of current branch segment are in second and third columns stored as doubles
			begin = REAL(eventSegs)[j + nrow * 1];
			end = REAL(eventSegs)[j + nrow * 2];
			
			//find next node to later check for shift point on branch
			if (j < nrow) {
				nnode = (int) REAL(eventSegs)[(j+1) + nrow * 0];
				nxtevent = (int) REAL(eventSegs)[(j+1) + nrow * 3];
			}
			
			//eventData is dataframe holding event parameters for the current tree
			//need to find the row that corresponds to the event index. in eventData
			//the rows are strictly ordered such that row 0 = event1, row 1 = event2, etc.
			Start = REAL(getListElement(eventData, "time"))[event-1];
			lam1 = REAL(getListElement(eventData, "lam1"))[event-1];
			lam2 = REAL(getListElement(eventData, "lam2"))[event-1];
			if (INTEGER(type)[0] == 0) {
			    mu1 = REAL(getListElement(eventData, "mu1"))[event-1];
			    mu2 = REAL(getListElement(eventData, "mu2"))[event-1];
			}
						
			//need to find which approximating segments match this branch segment
			//these are passed in strict order so we only need to search top to bottom
			//and can ignore everything we've been over already
			for (l = place_holder; l < nsegs; l++) {
				if ( (int) getDblMatrixELT(segmat, l, 0) == node) {
					//isGoodStart = REAL(segbegin)[l] >= begin;
					isGoodStart = ( (getDblMatrixELT(segmat, l, 1) - begin) >= 0. || ( (getDblMatrixELT(segmat, l, 1) - begin) < 0. && (getDblMatrixELT(segmat, l, 1) - begin) >= -1.*eps) );
					//isGoodEnd = REAL(segend)[l] <= end;
					isGoodEnd =  ( (getDblMatrixELT(segmat, l, 2) - end) <= 0. || ( (getDblMatrixELT(segmat, l, 2) - end) > 0. && (getDblMatrixELT(segmat, l, 2) - end) <= eps) );

					if (isGoodStart && isGoodEnd) {					
						relStart = getDblMatrixELT(segmat, l, 1) - Start;
						relEnd = getDblMatrixELT(segmat, l, 2) - Start;
						ret = getMeanRateExponential(relStart,relEnd,lam1,lam2);
						REAL(rates)[l] += ret/((double) nsamples);
						if (INTEGER(type)[0] == 0) {
							ret = getMeanRateExponential(relStart,relEnd,mu1,mu2);
							REAL(erates)[l] += ret/((double) nsamples);
						}
					}
					//check for shift straddle
					if (node == nnode) {
						isGoodStart = getDblMatrixELT(segmat, l, 1) < end;
						isGoodEnd = getDblMatrixELT(segmat, l, 2) > end;
						if (isGoodStart && isGoodEnd) {	
							relStart = getDblMatrixELT(segmat, l, 1) - Start;
							relEnd = end - Start;
							leftshift = getTimeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
							if (INTEGER(type)[0] == 0) {
								eleftshift = getTimeIntegratedBranchRate(relStart,relEnd,mu1,mu2);
							}
							relStart = 0.;
							relEnd = getDblMatrixELT(segmat, l, 2) - end;
							lam1 = REAL(getListElement(eventData, "lam1"))[nxtevent-1];
							lam2 = REAL(getListElement(eventData, "lam2"))[nxtevent-1];
							
							rightshift = getTimeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
							ret = (leftshift+rightshift)/(getDblMatrixELT(segmat, l, 2) - getDblMatrixELT(segmat, l, 1));
							REAL(rates)[l] += ret/((double) nsamples);
							
							if (INTEGER(type)[0] == 0) {
								mu1 = REAL(getListElement(eventData, "mu1"))[nxtevent-1];
								mu2 = REAL(getListElement(eventData, "mu2"))[nxtevent-1];
								erightshift = getTimeIntegratedBranchRate(relStart,relEnd,mu1,mu2);
								ret = (eleftshift+erightshift)/(getDblMatrixELT(segmat, l, 2) - getDblMatrixELT(segmat, l, 1));
								REAL(erates)[l] += ret/((double) nsamples);
							}
							place_holder = l; place_holder++;
							break;
						}
					}
				}
				else {
					place_holder = l;
					break;
				}
			}			
		}
		UNPROTECT(2); nprotect -= 2; //protected eventSegs and eventData, which we no longer need
	}
	if (INTEGER(type)[0] == 0) {
		SEXP retlist;
		PROTECT(retlist = allocVector(VECSXP, 2)); nprotect++;
		SET_VECTOR_ELT(retlist, 0, rates);
		SET_VECTOR_ELT(retlist, 1, erates);
		UNPROTECT(nprotect);
		return retlist;
	}
	UNPROTECT(nprotect);
	return rates;
}
