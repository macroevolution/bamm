\name{maximumShiftCredibilityTree}
\alias{maximumShiftCredibilityTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
maximumShiftCredibilityTree(ephy, maximize = "product")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{maximize}{
%%     ~~Describe \code{maximize} here~~
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
function (ephy, maximize = "product") 
{
    if (!"bammdata" \%in\% class(ephy)) {
        stop("Object ephy must be of class bammdata\n")
    }
    probvec <- numeric(length(ephy$eventData))
    mtree <- marginalShiftProbsTree(ephy)
    px <- mtree$edge.length
    for (i in 1:length(ephy$eventData)) {
        hasShift <- ephy$edge[, 2] \%in\% ephy$eventData[[i]]$node
        branchprobs <- (hasShift) * px + (!hasShift) * (1 - px)
        if (maximize == "product") {
            probvec[i] <- sum(log(branchprobs))
        }
        else if (maximize == "sum") {
            probvec[i] <- sum(branchprobs)
        }
        else {
            stop("Unsupported optimize criterion in maximumShiftCredibilityTree")
        }
    }
    best <- which(probvec == max(probvec))
    bestconfigs <- list()
    index <- 0
    while (length(best) > 0) {
        index <- index + 1
        lv <- logical(length = length(best))
        for (i in 1:length(best)) {
            lv[i] <- areEventConfigurationsIdentical(ephy, best[1], 
                best[i])
        }
        bestconfigs[[index]] <- best[lv]
        best <- best[!lv]
    }
    sampleindex <- numeric(length(bestconfigs))
    for (i in 1:length(bestconfigs)) {
        sampleindex[i] <- bestconfigs[[i]][1]
    }
    obj <- list()
    obj$bestconfigs <- bestconfigs
    obj$scores <- probvec
    obj$optimalityType = maximize
    obj$sampleindex <- sampleindex
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
