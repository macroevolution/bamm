\name{getBySpeciesRateMatrix}
\alias{getBySpeciesRateMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getBySpeciesRateMatrix(ephy, start.time, nbreaks, ndr = TRUE, node)
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
  \item{node}{
%%     ~~Describe \code{node} here~~
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
function (ephy, start.time, nbreaks, ndr = TRUE, node) 
{
    spset <- ephy$tip.label
    mm <- matrix(NA, nrow = length(spset), ncol = nbreaks)
    for (i in 1:nrow(mm)) {
        mm[i, ] <- getSpeciesRateThroughTime(ephy, start.time, 
            nbreaks, ndr, species = spset[i])
    }
    rownames(mm) <- spset
    return(mm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
