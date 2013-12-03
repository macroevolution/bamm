\name{segMap}
\alias{segMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
segMap(nodes, begin, end, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nodes}{
%%     ~~Describe \code{nodes} here~~
}
  \item{begin}{
%%     ~~Describe \code{begin} here~~
}
  \item{end}{
%%     ~~Describe \code{end} here~~
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
function (nodes, begin, end, tau) 
{
    foo = function(x, tau) {
        ret = seq(x[2], x[3], by = tau)
        if (length(ret) == 1) 
            return(matrix(x, nrow = 1))
        ret = seq(x[2], x[3], length.out = length(ret))
        ret = rep(ret, each = 2)
        ret = ret[2:(length(ret) - 1)]
        ret = matrix(ret, ncol = 2, byrow = TRUE)
        return(cbind(matrix(rep(as.integer(x[1]), nrow(ret)), 
            ncol = 1), ret))
    }
    times = cbind(nodes, begin, end)
    ret = apply(times, 1, foo, tau)
    return(do.call(rbind, ret))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
