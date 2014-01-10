\name{getmrca}
\alias{getmrca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find most recent common ancestors
}
\description{
Calculates the most recent common ancestor for each pair of tips.
}
\usage{
getmrca(phy, t1, t2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
an object of class \code{phylo}.
}
  \item{t1}{
a vector of mode integer or character corresponding to tips in \code{phy}.
}
  \item{t2}{
a vector of mode integer or character corresponding to tips in \code{phy}.
}
}
\details{
Finds the most recent common ancestor for each pair of tips where pairs are
defined as (\code{t1}[1], \code{t2}[1]), (\code{t1}[2], \code{t2}[2]), ... , (\code{t1}[i], \code{t2}[i]), ... ,(\code{t1}[n], \code{t2}[n]).
}
\value{
A vector of node numbers of the common ancestor for each pair of tips.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mike Grundler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
