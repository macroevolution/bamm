\name{getRateThroughTimeMatrix}
\alias{getRateThroughTimeMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getRateThroughTimeMatrix(ephy, start.time = NULL, end.time = NULL, nslices = 100, node = NULL, nodetype = "include")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{start.time}{
%%     ~~Describe \code{start.time} here~~
}
  \item{end.time}{
%%     ~~Describe \code{end.time} here~~
}
  \item{nslices}{
%%     ~~Describe \code{nslices} here~~
}
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{nodetype}{
%%     ~~Describe \code{nodetype} here~~
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
function (ephy, start.time = NULL, end.time = NULL, nslices = 100, 
    node = NULL, nodetype = "include") 
{
    if (!"bamm-data" \%in\% class(ephy)) {
        stop("Object ephy must be of class bamm-data\n")
    }
    if (is.null(node)) {
        nodeset <- ephy$edge[, 2]
    }
    else if (!is.null(node) & nodetype == "include") {
        nodeset <- getDesc(ephy, node)$desc_set
    }
    else if (!is.null(node) & nodetype == "exclude") {
        nodeset <- setdiff(ephy$edge[, 2], getDesc(ephy, node)$desc_set)
    }
    else {
        stop("error in getRateTHroughTimeMatrix\n")
    }
    bt <- branching.times(ephy)
    maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[, 
        1]))])
    if (is.null(start.time)) {
        start.time <- max(bt) - maxpossible
    }
    if (is.null(end.time)) {
        end.time <- max(bt)
    }
    tvec <- seq(start.time, end.time, length.out = nslices)
    mm <- matrix(NA, nrow = length(ephy$eventBranchSegs), ncol = length(tvec))
    mumat <- matrix(NA, nrow = length(ephy$eventBranchSegs), 
        ncol = length(tvec))
    for (i in 1:nrow(mm)) {
        es <- ephy$eventBranchSegs[[i]]
        events <- ephy$eventData[[i]]
        for (k in 1:length(tvec)) {
            isGoodTime <- es[, 2] <= tvec[k] & es[, 3] >= tvec[k]
            if (is.null(node)) {
                isGoodNode <- rep(TRUE, nrow(es))
            }
            else {
                isGoodNode <- es[, 1] \%in\% nodeset
            }
            estemp <- es[isGoodTime & isGoodNode, ]
            tvv <- tvec[k] - events$time[estemp[, 4]]
            rates <- exponentialRate(tvv, events$lam1[estemp[, 
                4]], events$lam2[estemp[, 4]])
            mm[i, k] <- mean(rates)
            mumat[i, k] <- mean(events$mu1[estemp[, 4]])
        }
    }
    obj <- list()
    obj$lambda <- mm
    obj$mu <- mumat
    obj$times <- tvec
    class(obj) <- "bamm-ratematrix"
    if (ephy$type == "diversification") {
        obj$type = "diversification"
    }
    else {
        obj$type = "traits"
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
