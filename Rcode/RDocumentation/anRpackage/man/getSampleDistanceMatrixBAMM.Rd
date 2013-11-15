\name{getSampleDistanceMatrixBAMM}
\alias{getSampleDistanceMatrixBAMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSampleDistanceMatrixBAMM(phylist, modeltree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylist}{
%%     ~~Describe \code{phylist} here~~
}
  \item{modeltree}{
%%     ~~Describe \code{modeltree} here~~
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
function (phylist, modeltree) 
{
    for (i in 1:length(phylist)) {
        phylist[[i]]$edge.length <- phylist[[i]]$edge.length * 
            modeltree$edge.length
    }
    dmat <- matrix(0, nrow = length(phylist), ncol = length(phylist))
    nsamples <- length(phylist)
    for (i in 1:(nsamples - 1)) {
        for (j in i:nsamples) {
            dmat[i, j] <- sum((phylist[[i]]$edge.length - phylist[[j]]$edge.length)^2)
        }
    }
    dmat <- dmat + t(dmat)
    return(dmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
