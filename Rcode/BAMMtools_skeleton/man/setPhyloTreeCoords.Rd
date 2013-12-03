\name{setPhyloTreeCoords}
\alias{setPhyloTreeCoords}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setPhyloTreeCoords(phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
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
function (phy) 
{
    phy = getStartStopTimes(phy)
    tH = as.numeric(max(branching.times(phy)))
    rootnd = as.integer(phy$Nnode + 2)
    ntip = as.integer(phy$Nnode + 1)
    anc = as.integer(phy$edge[, 1])
    desc = as.integer(phy$edge[, 2])
    nnode = as.integer(dim(phy$edge)[1] + 1)
    bl = as.numeric(phy$edge.length/tH)
    ndorder = .C("postorder_tree_traverse", anc, desc, rootnd, 
        nnode, integer(nnode))[[5]]
    ndorder = as.integer(ndorder)
    L = .C("setphylotreecoords", anc, desc, ndorder, bl, ntip, 
        rootnd, nnode, numeric(nrow(phy$edge) * 4), numeric(nrow(phy$edge) * 
            4), numeric(4))
    root = matrix(L[[10]], nrow = 1)
    xy = matrix(L[[9]], nrow = nrow(phy$edge), ncol = 4)
    bar = matrix(L[[8]], nrow = nrow(phy$edge), ncol = 4)
    xy = rbind(c(xy[1, 1], sum(xy[1:2, 4])/2, xy[1, 1], sum(xy[1:2, 
        4])/2), xy)
    bar = rbind(root, bar)
    rownames(xy) = c(phy$edge[1, 1], phy$edge[, 2])
    colnames(xy) = c("x0", "y0", "x1", "y1")
    return(list(segs = xy, arcs = bar))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
