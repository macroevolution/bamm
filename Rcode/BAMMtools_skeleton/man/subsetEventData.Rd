\name{subsetEventData}
\alias{subsetEventData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subsetEventData(ephy, index = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
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
function (ephy, index = NULL) 
{
    if (!"bamm-data" \%in\% class(ephy)) {
        stop("Object ephy must be of class bamm-data\n")
    }
    nsamples <- length(ephy$eventData)
    ss <- which((1:nsamples) \%in\% index)
    badsubset <- setdiff(index, 1:nsamples)
    if (length(badsubset) > 0) {
        cat("Bad call to BAMMtools::subsetEventData\n")
        cat("You have << ", nsamples, " >> samples in your bamm-data object \n")
        stop("Attempt to access invalid samples. Check index.")
    }
    obj <- list()
    obj$edge <- ephy$edge
    obj$Nnode <- ephy$Nnode
    obj$tip.label <- ephy$tip.label
    obj$edge.length <- ephy$edge.length
    obj$begin <- ephy$begin
    obj$end <- ephy$end
    obj$downseq <- ephy$downseq
    obj$lastvisit <- ephy$lastvisit
    obj$numberEvents <- ephy$numberEvents[ss]
    obj$eventData <- ephy$eventData[ss]
    obj$eventVectors <- ephy$eventVectors[ss]
    obj$tipStates <- ephy$tipStates[ss]
    obj$tipLambda <- ephy$tipLambda[ss]
    obj$tipMu <- ephy$tipMu[ss]
    obj$eventBranchSegs <- ephy$eventBranchSegs[ss]
    obj$meanTipLambda <- ephy$meanTipLambda
    obj$meanTipMu <- ephy$meanTipMu
    obj$type <- ephy$type
    class(obj) <- c("phylo", "bamm-data")
    attributes(obj)$order = attributes(ephy)$order
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
