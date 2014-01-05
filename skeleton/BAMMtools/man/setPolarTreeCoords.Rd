\name{setPolarTreeCoords}
\alias{setPolarTreeCoords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setPolarTreeCoords(phy, vtheta, rbf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
}
  \item{vtheta}{
%%     ~~Describe \code{vtheta} here~~
}
  \item{rbf}{
%%     ~~Describe \code{rbf} here~~
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
function (phy, vtheta, rbf) 
{
    phy = getStartStopTimes(phy)
    tH = max(branching.times(phy))
    rootnd = as.integer(phy$Nnode + 2)
    ntip = as.integer(phy$Nnode + 1)
    anc = as.integer(phy$edge[, 1])
    desc = as.integer(phy$edge[, 2])
    nnode = as.integer(dim(phy$edge)[1] + 1)
    vtheta = as.numeric(vtheta * (pi/180))
    ths = as.numeric((2 * pi - vtheta)/phy$Nnode)
    ndorder = .C("postorder_tree_traverse", anc, desc, rootnd, 
        nnode, integer(nnode))[[5]]
    ndorder = as.integer(ndorder)
    L = .C("setpolartreecoords", anc, desc, ndorder, ntip, rootnd, 
        nnode, ths, numeric(nrow(phy$edge) * 3), numeric(3))
    root = matrix(L[[9]], nrow = 1)
    theta = matrix(L[[8]], nrow = nrow(phy$edge), ncol = 3)
    theta = rbind(root, theta)
    rb = tH * rbf
    x0 = c(rb, rb + (phy$begin/tH)) * cos(theta[, 1])
    y0 = c(rb, rb + (phy$begin/tH)) * sin(theta[, 1])
    x1 = c(rb, rb + (phy$end/tH)) * cos(theta[, 1])
    y1 = c(rb, rb + (phy$end/tH)) * sin(theta[, 1])
    ret = cbind(x0, y0, x1, y1, theta[, 1])
    rownames(ret) = c(phy$edge[1, 1], phy$edge[, 2])
    colnames(ret) = c("x0", "y0", "x1", "y1", "theta")
    return(list(segs = ret, arcs = theta[, 2:3]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
