\name{getMeanTipRateStates}
\alias{getMeanTipRateStates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMeanTipRateStates(ephy, use.names = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{use.names}{
%%     ~~Describe \code{use.names} here~~
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
function (ephy, use.names = FALSE) 
{
    if (!"bamm-data" \%in\% class(ephy)) {
        stop("Object ephy must be of class bamm-data\n")
    }
    if (use.names == FALSE) {
        return(ephy$meanTipLambda - ephy$meanTipMu)
    }
    else {
        ret <- ephy$meanTipLambda - ephy$meanTipMu
        names(ret) <- ephy$tip.label
        return(ret)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
