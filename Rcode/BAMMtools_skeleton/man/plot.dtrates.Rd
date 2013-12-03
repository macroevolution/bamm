\name{plot.dtrates}
\alias{plot.dtrates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.dtrates(ephy, method = "phylogram", tau = 0.01, index = NULL, vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, hrates = FALSE, lwd = 3, cex = 1, ncolors = 64, pal = "Spectral", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{tau}{
%%     ~~Describe \code{tau} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{vtheta}{
%%     ~~Describe \code{vtheta} here~~
}
  \item{rbf}{
%%     ~~Describe \code{rbf} here~~
}
  \item{show}{
%%     ~~Describe \code{show} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{hrates}{
%%     ~~Describe \code{hrates} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{ncolors}{
%%     ~~Describe \code{ncolors} here~~
}
  \item{pal}{
%%     ~~Describe \code{pal} here~~
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
function (ephy, method = "phylogram", tau = 0.01, index = NULL, 
    vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, hrates = FALSE, 
    lwd = 3, cex = 1, ncolors = 64, pal = "Spectral", ...) 
{
    if ("bamm-data" \%in\% class(ephy)) 
        phy = as.phylo.bammdata(ephy)
    else stop("Trying to work with a non-bammdata object")
    if (!is.binary.tree(phy)) 
        stop("Function requires fully bifurcating tree.")
    if (any(phy$edge.length == 0)) 
        warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero")
    if (!("dtrates" \%in\% names(ephy))) 
        ephy = dtRates(ephy, tau, index)
    if (sum(is.na(ephy$dtrates$rates))) {
        warning(sprintf("Found \%s NA values. Coercing to zero.", 
            sum(is.na(ephy$dtrates$rates))))
        ephy$dtrates$rates[is.na(ephy$dtrates$rates)] = 0
    }
    if (method == "polar") {
        phy = getStartStopTimes(phy)
        ret = setPolarTreeCoords(phy, vtheta, rbf)
        tH = max(branching.times(phy))
        rb = tH * rbf
    }
    else if (method == "phylogram") {
        ret = setPhyloTreeCoords(phy)
    }
    else {
        stop("Unimplemented method")
    }
    x0 = ret$segs[, 1]
    y0 = ret$segs[, 2]
    x1 = ret$segs[, 3]
    y1 = ret$segs[, 4]
    tau = ephy$dtrates$tau
    edge.color = colorMap(ephy$dtrates$rates, pal, ncolors)
    p = cbind(x0[-1], y0[-1], x1[-1], y1[-1], phy$edge[, 2])
    p = apply(p, 1, mkdtsegs, tau)
    p = do.call(rbind, p)
    x0 = c(x0[1], p[, 1])
    x1 = c(x1[1], p[, 2])
    y0 = c(y0[1], p[, 3])
    y1 = c(y1[1], p[, 4])
    offset = table(p[, 5])[as.character(unique(p[, 5]))]
    arc.color = c(edge.color[1], edge.color[match(unique(p[, 
        5]), p[, 5]) + offset])
    edge.color = c(edge.color[1], edge.color)
    if (show) {
        op = par(no.readonly = TRUE)
        par(...)
        plot.new()
        ofs = 0
        if (labels) {
            ofs = max(nchar(phy$tip.label) * 0.03 * cex)
        }
        if (method == "polar") {
            plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                ofs), ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                ofs), asp = 1)
            segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, 
                lend = 2)
            arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + 
                phy$end/tH), border = arc.color, lwd = lwd)
            if (labels) {
                for (k in 1:length(phy$tip.label)) {
                  text(ret$segs[-1, ][phy$edge[, 2] == k, 3], 
                    ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k], 
                    cex = cex, srt = (180/pi) * ret$arcs[-1, 
                      ][phy$edge[, 2] == k, 1], adj = c(0, NA))
                }
            }
        }
        if (method == "phylogram") {
            plot.window(xlim = c(0, 1 + ofs), ylim = c(0, phy$Nnode * 
                1/(phy$Nnode + 1)), asp = 1)
            segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, 
                lend = 2)
            isTip = phy$edge[, 2] <= phy$Nnode + 1
            isTip = c(FALSE, isTip)
            segments(ret$arcs[!isTip, 1], ret$arcs[!isTip, 2], 
                ret$arcs[!isTip, 3], ret$arcs[!isTip, 4], col = arc.color[!isTip], 
                lwd = lwd, lend = 2)
            if (labels) {
                text(ret$segs[-1, ][phy$edge[, 2] <= phy$Nnode + 
                  1, 3], ret$segs[-1, ][phy$edge[, 2] <= phy$Nnode + 
                  1, 4], phy$tip.label, cex = cex, pos = 4, offset = 0.25)
            }
        }
        if (hrates) {
            histRates(ephy$dtrates$rates, pal, ncolors)
        }
    }
    index = order(as.numeric(rownames(ret$segs)))
    if (method == "phylogram") {
        assign("last_plot.phylo", list(Ntip = phy$Nnode + 1, 
            Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                3], yy = ret$segs[index, 4]), envir = .PlotPhyloEnv)
    }
    else if (method == "polar") {
        assign("last_plot.phylo", list(Ntip = phy$Nnode + 1, 
            Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                3], yy = ret$segs[index, 4], theta = ret$segs[index, 
                5], rb = rb), envir = .PlotPhyloEnv)
    }
    par(op)
    invisible(ret$segs[-1, ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
