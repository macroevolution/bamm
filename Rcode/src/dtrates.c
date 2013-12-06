#include <R.h>
#include <Rinternals.h>

//forward function declarations
SEXP getListElement(SEXP list, char *str);
SEXP dtrates(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample);
//SEXP getMatrixColumn(SEXP matrix, int col);
double getDblMatrixELT(SEXP matrix, int row, int col);
double getMeanRateExponential(double t1, double t2, double p1, double p2);
double getTimeIntegratedBranchRate(double t1, double t2, double p1, double p2);

SEXP getListElement(SEXP list, char *str)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < LENGTH(list); i++)
	{
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
		{
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}

/***
SEXP getMatrixColumn(SEXP matrix, int col)
{
	int nrow = INTEGER(getAttrib(matrix, R_DimSymbol))[0];
	
	SEXP ret;
	ret = PROTECT(allocVector(REALSXP, nrow));
	
	int i;
	for (i = 0; i < nrow; i++)
	{
		REAL(ret)[i] = REAL(matrix)[i + nrow * col];
	}
	return ret;
}
***/

double getDblMatrixELT(SEXP matrix, int row, int col)
{
	int nrow = INTEGER(getAttrib(matrix, R_DimSymbol))[0];
	return REAL(matrix)[row + nrow*col];
}

double getMeanRateExponential(double t1, double t2, double p1, double p2)
{	
	if (p2 == 0.) 
	{
		return p1;
	}
	return (p1/p2)*(exp(t2*p2) - exp(t1*p2))/(t2 - t1);
}

double getTimeIntegratedBranchRate(double t1, double t2, double p1, double p2)
{
	if (p2 == 0.) 
	{
		return (t2 - t1) * p1;
	}
	return (p1/p2)*(exp(t2*p2) - exp(t1*p2));
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

#define segmat(l,k) (getDblMatrixELT(segmat, l, k))

SEXP dtrates(SEXP ephy, SEXP segmat, SEXP tol, SEXP sample)
{
	double eps = REAL(tol)[0];
	
	int k, j, l, nprotect = 0;
	//int nsamples = LENGTH(getListElement(ephy, "eventBranchSegs"));
	int nsamples = LENGTH(sample);
	
	/***
	SEXP nodeseg = getMatrixColumn(segmat, 0); nprotect++;	
	SEXP segbegin = getMatrixColumn(segmat, 1); nprotect++;
	SEXP segend = getMatrixColumn(segmat, 2); nprotect++;
	int nsegs = LENGTH(nodeseg);
	***/
	
	int nsegs = INTEGER(getAttrib(segmat, R_DimSymbol))[0];
	
	SEXP rates;
	PROTECT(rates = allocVector(REALSXP, nsegs)); nprotect++;
	
	for(k = 0; k < nsegs; k++)
	{
		REAL(rates)[k] = 0.;
	}
	
	int nrow, node, nnode, event, nxtevent, isGoodStart, isGoodEnd, place_holder;
	double begin, end, Start, lam1, lam2, relStart, relEnd, rightshift, leftshift, ret;		
	//for (k = 0; k < nsamples; k++)
	for (k = INTEGER(sample)[0] - 1; k < INTEGER(sample)[nsamples - 1]; k++)
	{
		SEXP eventSegs, eventData;
		
		eventSegs = PROTECT(VECTOR_ELT(getListElement(ephy, "eventBranchSegs"), k)); nprotect++;
		eventData = PROTECT(VECTOR_ELT(getListElement(ephy, "eventData"), k)); nprotect++;
				
		nrow = INTEGER(getAttrib(eventSegs, R_DimSymbol))[0];
		place_holder = 0;
		for(j = 0; j < nrow; j++) //move down the rows of eventSegs
		{
			//eventSegs is 4 column matrix, node is in first column stored as double
			node = (int) REAL(eventSegs)[j + nrow * 0];
			
			//event index is in fourth column stored as double
			event = (int) REAL(eventSegs)[j  + nrow * 3];
			
			//begin and end of current branch segment are in second and third columns stored as doubles
			begin = REAL(eventSegs)[j + nrow * 1];
			end = REAL(eventSegs)[j + nrow * 2];
			
			//find next node to later check for shift point on branch
			if (j < nrow)
			{
				nnode = (int) REAL(eventSegs)[(j+1) + nrow * 0];
				nxtevent = (int) REAL(eventSegs)[(j+1) + nrow * 3];
			}
			
			//eventData is dataframe holding event parameters for the current tree
			//need to find the row that corresponds to the event index. in eventData
			//the rows are strictly ordered such that row 0 = event1, row 1 = event2, etc.
			Start = REAL(getListElement(eventData, "time"))[event-1];
			lam1 = REAL(getListElement(eventData, "lam1"))[event-1];
			lam2 = REAL(getListElement(eventData, "lam2"))[event-1];
						
			//need to find which approximating segments match this branch segment
			//these are passed in strict order so we only need to search top to bottom
			//and can ignore everything we've been over already
			for (l = place_holder; l < nsegs; l++)
			{
				if ( (int) segmat(l, 0) == node)
				{
					//isGoodStart = REAL(segbegin)[l] >= begin;
					isGoodStart = ( (segmat(l,1) - begin) >= 0. || ( (segmat(l,1) - begin) < 0. && (segmat(l,1) - begin) >= -1.*eps) );
					//isGoodEnd = REAL(segend)[l] <= end;
					isGoodEnd =  ( (segmat(l,2) - end) <= 0. || ( (segmat(l,2) - end) > 0. && (segmat(l,2) - end) <= eps) );

					if (isGoodStart && isGoodEnd)
					{					
						relStart = segmat(l,1) - Start;
						relEnd = segmat(l,2) - Start;
						ret = getMeanRateExponential(relStart,relEnd,lam1,lam2);
						
						REAL(rates)[l] += ret/((double) nsamples);
					}
					//check for shift straddle
					if (node == nnode)
					{
						isGoodStart = segmat(l,1) < end;
						isGoodEnd = segmat(l,2) > end;
						if (isGoodStart && isGoodEnd)
						{	
							relStart = segmat(l,1) - Start;
							relEnd = end - Start;
							leftshift = getTimeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
							
							relStart = 0.;
							relEnd = segmat(l,2) - end;
							lam1 = REAL(getListElement(eventData, "lam1"))[nxtevent-1];
							lam2 = REAL(getListElement(eventData, "lam2"))[nxtevent-1];
							rightshift = getTimeIntegratedBranchRate(relStart,relEnd,lam1,lam2);
							
							ret = (leftshift+rightshift)/(segmat(l,2) - segmat(l,1));
							
							REAL(rates)[l] += ret/((double) nsamples);
							place_holder = l; place_holder++;
							break;
						}
					}
				}
				else
				{
					place_holder = l;
					break;
				}
			}			
		}
		UNPROTECT(2); nprotect -= 2; //protected eventSegs and eventData, which we no longer need
	}
	UNPROTECT(nprotect);
	return rates;
}