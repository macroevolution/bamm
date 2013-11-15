\name{getCladeRates}
\alias{getCladeRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCladeRates(ephy, weights = "branchlengths", node = NULL, nodetype = "include", verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{nodetype}{
%%     ~~Describe \code{nodetype} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (ephy, weights = "branchlengths", node = NULL, nodetype = "include", 
    verbose = F) 
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
    lambda_vector <- numeric(length(ephy$eventBranchSegs))
    mu_vector <- numeric(length(ephy$eventBranchSegs))
    for (i in 1:length(ephy$eventBranchSegs)) {
        if (verbose) {
            cat("Processing sample ", i, "\n")
        }
        esegs <- ephy$eventBranchSegs[[i]]
        esegs <- esegs[esegs[, 1] \%in\% nodeset, ]
        events <- ephy$eventData[[i]]
        events <- events[order(events$index), ]
        relsegmentstart <- esegs[, 2] - ephy$eventData[[i]]$time[esegs[, 
            4]]
        relsegmentend <- esegs[, 3] - ephy$eventData[[i]]$time[esegs[, 
            4]]
        lam1 <- ephy$eventData[[i]]$lam1[esegs[, 4]]
        lam2 <- ephy$eventData[[i]]$lam2[esegs[, 4]]
        mu1 <- ephy$eventData[[i]]$mu1[esegs[, 4]]
        mu2 <- ephy$eventData[[i]]$mu2[esegs[, 4]]
        seglengths <- esegs[, 3] - esegs[, 2]
        wts <- seglengths/sum(seglengths)
        lamseg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, 
            lam1, lam2)/seglengths
        museg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, 
            mu1, mu2)/seglengths
        lambda_vector[i] <- sum(lamseg * wts)
        mu_vector[i] <- sum(museg * wts)
    }
    return(list(lambda = lambda_vector, mu = mu_vector))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
