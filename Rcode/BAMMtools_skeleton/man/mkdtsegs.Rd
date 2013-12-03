\name{mkdtsegs}
\alias{mkdtsegs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mkdtsegs(x, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
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
function (x, tau) 
{
    bn = sqrt((x[3] - x[1])^2 + (x[4] - x[2])^2)
    len = bn/tau
    if (len\%\%1 == 0) 
        len = len + 1
    j = seq(x[1], x[3], length.out = len)
    if (length(j) == 1) 
        return(matrix(x[c(1, 3, 2, 4, 5)], nrow = 1))
    k = seq(x[2], x[4], length.out = len)
    j = rep(j, each = 2)
    j = j[2:(length(j) - 1)]
    j = matrix(j, ncol = 2, byrow = TRUE)
    k = rep(k, each = 2)
    k = k[2:(length(k) - 1)]
    k = matrix(k, ncol = 2, byrow = TRUE)
    l = matrix(rep(x[5], nrow(j)), ncol = 1)
    return(cbind(j, k, l))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
