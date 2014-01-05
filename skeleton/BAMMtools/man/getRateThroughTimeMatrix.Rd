\name{getRateThroughTimeMatrix}
\alias{getRateThroughTimeMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp
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
    if (!"bammdata" \%in\% class(ephy)) {
        stop("Object ephy must be of class 'bammdata'\n")
    }
    if (is.null(node)) {
        nodeset <- ephy$edge[, 2]
    }
    else if (!is.null(node) & nodetype == "include") {
        nodeset <- unlist(sapply(node, function(x) getDesc(ephy, 
            x)$desc_set))
    }
    else if (!is.null(node) & nodetype == "exclude") {
        nodeset <- setdiff(ephy$edge[, 2], unlist(sapply(node, 
            function(x) getDesc(ephy, x)$desc_set)))
    }
    else {
        stop("error in getRateThroughTimeMatrix\n")
    }
    bt <- branching.times(as.phylo.bammdata(ephy))
    maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[, 
        1]))])
    if (is.null(start.time)) {
        start.time <- max(bt) - maxpossible
    }
    if (is.null(end.time)) {
        end.time <- max(bt)
    }
    tvec <- seq(start.time, end.time, length.out = nslices)
    tol = 1 * 10^-decimals(ephy$eventBranchSegs[[1]][1, 2])
    mm <- matrix(NA, nrow = length(ephy$eventBranchSegs), ncol = length(tvec))
    mumat <- matrix(NA, nrow = length(ephy$eventBranchSegs), 
        ncol = length(tvec))
    goodTime <- function(vec, val, tol) {
        (vec[, 2] - val <= tol) & (val - vec[, 3] <= tol)
    }
    for (i in 1:nrow(mm)) {
        es <- ephy$eventBranchSegs[[i]]
        events <- ephy$eventData[[i]]
        isGoodNode <- rep(TRUE, nrow(es))
        for (k in 1:length(tvec)) {
            isGoodTime <- goodTime(es, tvec[k], tol = tol)
            if (!(is.null(node))) {
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
    if (ephy$type == "diversification") {
        obj$lambda <- mm
        obj$mu <- mumat
    }
    if (ephy$type == "trait") {
        obj$beta <- mm
    }
    obj$times <- tvec
    class(obj) <- "bamm-ratematrix"
    if (ephy$type == "diversification") {
        obj$type = "diversification"
    }
    else {
        obj$type = "trait"
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
