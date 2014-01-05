\name{addBAMMshifts}
\alias{addBAMMshifts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add BAMM shifts to a plotted tree
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addBAMMshifts(ephy, method, index, cex = 1, pch = 21, col = 1, bg = 2, multi = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{bg}{
%%     ~~Describe \code{bg} here~~
}
  \item{multi}{
%%     ~~Describe \code{multi} here~~
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
function (ephy, method, index, cex = 1, pch = 21, col = 1, bg = 2, 
    multi = FALSE) 
{
    if (!"bammdata" \%in\% class(ephy)) 
        stop("Object ephy must be of class bammdata")
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (!multi) {
        op = par(no.readonly = TRUE)
        par(lastPP$pp)
    }
    shiftnodes = getShiftNodesFromIndex(ephy, index)
    isShift = ephy$eventData[[index]]$node \%in\% shiftnodes
    times = ephy$eventData[[index]]$time[isShift]
    if (method == "phylogram") {
        if (max(lastPP$xx) <= 1) {
            XX = times/max(branching.times(as.phylo.bammdata(ephy)))
        }
        else {
            XX = times
        }
        YY = lastPP$yy[shiftnodes]
    }
    else if (method == "polar") {
        rb = lastPP$rb
        XX = (rb + times/max(branching.times(as.phylo.bammdata(ephy)))) * 
            cos(lastPP$theta[shiftnodes])
        YY = (rb + times/max(branching.times(as.phylo.bammdata(ephy)))) * 
            sin(lastPP$theta[shiftnodes])
    }
    points(XX, YY, pch = pch, cex = cex, col = col, bg = bg)
    if (!multi) {
        par(op)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
