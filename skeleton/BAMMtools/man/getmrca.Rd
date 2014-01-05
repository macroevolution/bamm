\name{getmrca}
\alias{getmrca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getmrca(phy, t1, t2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
}
  \item{t1}{
%%     ~~Describe \code{t1} here~~
}
  \item{t2}{
%%     ~~Describe \code{t2} here~~
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
function (phy, t1, t2) 
{
    ne = as.integer(dim(phy$edge)[1])
    npair = as.integer(length(t1))
    anc = as.integer(phy$edge[, 1])
    desc = as.integer(phy$edge[, 2])
    root = as.integer(length(phy$tip.label) + 1)
    .C("fetchmrca", anc, desc, root, ne, npair, as.integer(t1), 
        as.integer(t2), integer(npair))[[8]]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
