\name{getMarginalBranchRateMatrix}
\alias{getMarginalBranchRateMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMarginalBranchRateMatrix(obj, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
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
function (obj, verbose = FALSE) 
{
    if (!"bamm-data" \%in\% class(obj) & class(obj[[1]]) != "phylo") {
        stop("Object must be of class bamm-data or a list of trees.\n")
    }
    if ("bamm-data" \%in\% class(obj)) {
        lammat <- matrix(0, ncol = length(obj$eventBranchSegs), 
            nrow = nrow(obj$edge))
        mumat <- matrix(0, ncol = length(obj$eventBranchSegs), 
            nrow = nrow(obj$edge))
        for (i in 1:length(obj$eventBranchSegs)) {
            if (verbose) {
                cat("Processing sample ", i, "\n")
            }
            esegs <- obj$eventBranchSegs[[i]]
            events <- obj$eventData[[i]]
            events <- events[order(events$index), ]
            relsegmentstart <- esegs[, 2] - obj$eventData[[i]]$time[esegs[, 
                4]]
            relsegmentend <- esegs[, 3] - obj$eventData[[i]]$time[esegs[, 
                4]]
            lam1 <- obj$eventData[[i]]$lam1[esegs[, 4]]
            lam2 <- obj$eventData[[i]]$lam2[esegs[, 4]]
            mu1 <- obj$eventData[[i]]$mu1[esegs[, 4]]
            mu2 <- obj$eventData[[i]]$mu2[esegs[, 4]]
            lamint <- timeIntegratedBranchRate(relsegmentstart, 
                relsegmentend, lam1, lam2)
            muint <- timeIntegratedBranchRate(relsegmentstart, 
                relsegmentend, mu1, mu2)
            seglengths <- esegs[, 3] - esegs[, 2]
            for (k in 1:nrow(obj$edge)) {
                isRightBranch <- esegs[, 1] == obj$edge[k, 2]
                lammat[k, i] <- sum(lamint[isRightBranch])/sum(seglengths[isRightBranch])
                mumat[k, i] <- sum(muint[isRightBranch])/sum(seglengths[isRightBranch])
            }
        }
        if (obj$type == "diversification") {
            return(list(lambda_branch_matrix = lammat, mu_branch_matrix = mumat))
        }
        if (obj$type == "traits") {
            return(list(beta_branch_matrix = lammat))
        }
    }
    else if (class(obj[[1]]) == "phylo") {
        return(sapply(obj, with, edge.length))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
