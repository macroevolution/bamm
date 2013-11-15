\name{getSpeciesRateThroughTime}
\alias{getSpeciesRateThroughTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSpeciesRateThroughTime(ephy, start.time, nbreaks = 10, ndr = TRUE, species)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{start.time}{
%%     ~~Describe \code{start.time} here~~
}
  \item{nbreaks}{
%%     ~~Describe \code{nbreaks} here~~
}
  \item{ndr}{
%%     ~~Describe \code{ndr} here~~
}
  \item{species}{
%%     ~~Describe \code{species} here~~
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
function (ephy, start.time, nbreaks = 10, ndr = TRUE, species) 
{
    if (!"bamm-data" \%in\% class(ephy)) {
        stop("Object ephy must be of class bamm-data\n")
    }
    tend <- max(branching.times(ephy)) * 0.999
    tstart <- start.time
    tseq <- seq(tstart, tend, length.out = nbreaks)
    res <- numeric(nbreaks)
    index <- which(ephy$tip.label == species)
    path <- getPathToRoot(ephy, node = index)
    for (k in 1:length(ephy$eventBranchSegs)) {
        ed <- ephy$eventData[[k]]
        for (z in 1:length(tseq)) {
            isGoodNode <- ephy$eventBranchSegs[[k]][, 1] \%in\% 
                path
            isGoodStart <- ephy$eventBranchSegs[[k]][, 2] <= 
                tseq[z]
            isGoodEnd <- ephy$eventBranchSegs[[k]][, 3] >= tseq[z]
            ev <- ephy$eventBranchSegs[[k]][, 4][isGoodNode & 
                isGoodStart & isGoodEnd]
            if (ndr) {
                lam <- exponentialRate(tseq[z] - ed$time[ev], 
                  ed$lam1[ev], ed$lam2[ev])
                mu <- ed$mu1[ev]
                res[z] <- res[z] + (lam - mu)
            }
            else {
                res[z] <- res[z] + exponentialRate(tseq[z] - 
                  ed$time[ev], ed$lam1[ev], ed$lam2[ev])
            }
        }
    }
    res <- res/length(ephy$eventBranchSegs)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
