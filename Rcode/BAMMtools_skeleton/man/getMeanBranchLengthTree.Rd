\name{getMeanBranchLengthTree}
\alias{getMeanBranchLengthTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMeanBranchLengthTree(obj, ndr = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
%%     ~~Describe \code{obj} here~~
}
  \item{ndr}{
%%     ~~Describe \code{ndr} here~~
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
function (obj, ndr = TRUE) 
{
    if ("bamm-data" \%in\% class(obj)) {
        v <- list(edge = obj$edge, Nnode = obj$Nnode, tip.label = obj$tip.label, 
            edge.length = obj$edge.length)
        attributes(v) <- list(names = c("edge", "Nnode", "tip.label", 
            "edge.length"), class = "phylo", order = "cladewise")
        obj <- getMarginalBranchRateMatrix(obj, verbose = FALSE)
        if (ndr) {
            el <- rowMeans(obj$lambda) - rowMeans(obj$mu)
        }
        else {
            el <- rowMeans(obj$lambda)
        }
        v$edge.length <- el
        obj <- list()
        obj$phy <- v
        obj$median <- median(el)
        obj$mean <- mean(el)
    }
    else if (class(obj[[1]]) == "phylo") {
        el = rowMeans(sapply(obj, with, edge.length))
        obj <- list()
        obj$phy <- obj[[1]]
        obj$median <- median(el)
        obj$mean <- mean(el)
    }
    else {
        stop("Method not implemented for supplied object class")
    }
    return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
