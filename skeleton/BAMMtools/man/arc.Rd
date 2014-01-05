\name{arc}
\alias{arc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Temporary title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
arc(x, y, theta1, theta2, rad, border, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{theta1}{
%%     ~~Describe \code{theta1} here~~
}
  \item{theta2}{
%%     ~~Describe \code{theta2} here~~
}
  \item{rad}{
%%     ~~Describe \code{rad} here~~
}
  \item{border}{
%%     ~~Describe \code{border} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, y, theta1, theta2, rad, border, ...) 
{
    steps = (theta2 - theta1)/100
    noTips = steps > 0
    steps = steps[noTips]
    theta1 = theta1[noTips]
    theta2 = theta2[noTips]
    rad = rad[noTips]
    border = border[noTips]
    for (i in 1:length(steps)) {
        xv = x + rad[i] * cos(seq(theta1[i], theta2[i], steps[i]))
        yv = y + rad[i] * sin(seq(theta1[i], theta2[i], steps[i]))
        lines(xv, yv, lend = 2, col = border[i], ...)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
