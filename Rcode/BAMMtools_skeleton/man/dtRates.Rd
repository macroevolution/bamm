\name{dtRates}
\alias{dtRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dtRates(ephy, tau, ism = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{ism}{
%%     ~~Describe \code{ism} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ephy, tau, ism = NULL) 
{
    if (!"bamm-data" \%in\% class(ephy)) {
        stop("Function requires a bammdata object")
    }
    ephy$eventBranchSegs = lapply(ephy$eventBranchSegs, function(x) x[order(x[, 
        1]), ])
    phy = as.phylo.bammdata(ephy)
    phy = getStartStopTimes(phy)
    if (attributes(phy)$order != "cladewise") {
        phy = reorder(phy, "cladewise")
    }
    tH = max(branching.times(phy))
    segmat = segMap(phy$edge[, 2], phy$begin/tH, phy$end/tH, 
        tau)
    segmat[, 2] = segmat[, 2] * tH
    segmat[, 3] = segmat[, 3] * tH
    tol = 1 * 10^-decimals(ephy$eventBranchSegs[[1]][1, 2])
    if (storage.mode(segmat) != "double") 
        stop("Exiting")
    if (storage.mode(tol) != "double") 
        stop("Exiting")
    if (storage.mode(ephy) != "list") 
        stop("Exiting")
    if (is.null(ism)) 
        ism = as.integer(1:length(ephy$eventBranchSegs))
    else ism = as.integer(ism)
    if (ism[length(ism)] > length(ephy$eventBranchSegs)) {
        warning("Sample index out of range")
        ism = as.integer(1:length(ephy$eventBranchSegs))
    }
    index = 1:nrow(segmat)
    rownames(segmat) = index
    segmat = segmat[order(segmat[, 1]), ]
    dtrates = .Call("dtrates", ephy, segmat, tol, ism)
    names(dtrates) = rownames(segmat)
    dtrates = dtrates[as.character(index)]
    names(dtrates) = NULL
    ephy$dtrates = list(tau = tau, rates = dtrates)
    return(ephy)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
