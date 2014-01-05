\name{safeCompare.old}
\alias{safeCompare.old}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
safeCompare.old(vec, val, FUN, tol = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
%%     ~~Describe \code{vec} here~~
}
  \item{val}{
%%     ~~Describe \code{val} here~~
}
  \item{FUN}{
%%     ~~Describe \code{FUN} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
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
function (vec, val, FUN, tol = 1e-04) 
{
    if (FUN == ">=") {
        ret = rep(FALSE, length(vec))
        ret[(vec - val) >= 0] = TRUE
        ret[abs(vec - val) <= tol] = TRUE
        return(ret)
    }
    if (FUN == "<=") {
        ret = rep(FALSE, length(vec))
        ret[(vec - val) <= 0] = TRUE
        ret[(vec - val) <= tol] = TRUE
        return(ret)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
