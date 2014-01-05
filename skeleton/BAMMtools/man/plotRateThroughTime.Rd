\name{plotRateThroughTime}
\alias{plotRateThroughTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotRateThroughTime(ephy, useMedian = FALSE, intervals = seq(from = 0, to = 1, by = 0.01), ratetype = "auto", nBins = 100, smooth = FALSE, smoothParam = 0.2, opacity = 0.01, intervalCol = "blue", avgCol = "red", start.time = NULL, end.time = NULL, node = NULL, nodetype = "include", plot = TRUE, cex.axis = 1, cex = 1.3, xline = 3.5, yline = 3.5, mar = c(6, 6, 1, 1), xticks = 5, yticks = 5, xlim = "auto", ylim = "auto", add = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{useMedian}{
%%     ~~Describe \code{useMedian} here~~
}
  \item{intervals}{
%%     ~~Describe \code{intervals} here~~
}
  \item{ratetype}{
%%     ~~Describe \code{ratetype} here~~
}
  \item{nBins}{
%%     ~~Describe \code{nBins} here~~
}
  \item{smooth}{
%%     ~~Describe \code{smooth} here~~
}
  \item{smoothParam}{
%%     ~~Describe \code{smoothParam} here~~
}
  \item{opacity}{
%%     ~~Describe \code{opacity} here~~
}
  \item{intervalCol}{
%%     ~~Describe \code{intervalCol} here~~
}
  \item{avgCol}{
%%     ~~Describe \code{avgCol} here~~
}
  \item{start.time}{
%%     ~~Describe \code{start.time} here~~
}
  \item{end.time}{
%%     ~~Describe \code{end.time} here~~
}
  \item{node}{
%%     ~~Describe \code{node} here~~
}
  \item{nodetype}{
%%     ~~Describe \code{nodetype} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{cex.axis}{
%%     ~~Describe \code{cex.axis} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{xline}{
%%     ~~Describe \code{xline} here~~
}
  \item{yline}{
%%     ~~Describe \code{yline} here~~
}
  \item{mar}{
%%     ~~Describe \code{mar} here~~
}
  \item{xticks}{
%%     ~~Describe \code{xticks} here~~
}
  \item{yticks}{
%%     ~~Describe \code{yticks} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{add}{
%%     ~~Describe \code{add} here~~
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
function (ephy, useMedian = FALSE, intervals = seq(from = 0, 
    to = 1, by = 0.01), ratetype = "auto", nBins = 100, smooth = FALSE, 
    smoothParam = 0.2, opacity = 0.01, intervalCol = "blue", 
    avgCol = "red", start.time = NULL, end.time = NULL, node = NULL, 
    nodetype = "include", plot = TRUE, cex.axis = 1, cex = 1.3, 
    xline = 3.5, yline = 3.5, mar = c(6, 6, 1, 1), xticks = 5, 
    yticks = 5, xlim = "auto", ylim = "auto", add = FALSE) 
{
    if (!any(c("bammdata", "bamm-ratematrix") \%in\% class(ephy))) {
        stop("ERROR: Object ephy must be of class 'bammdata' or 'bamm-ratematrix'.\n")
    }
    if (!is.logical(useMedian)) {
        stop("ERROR: useMedian must be either TRUE or FALSE.")
    }
    if (!any(c("numeric", "NULL") \%in\% class(intervals))) {
        stop("ERROR: intervals must be either 'NULL' or a vector of quantiles.")
    }
    if (!is.logical(smooth)) {
        stop("ERROR: smooth must be either TRUE or FALSE.")
    }
    if ("bammdata" \%in\% class(ephy)) {
        rmat <- getRateThroughTimeMatrix(ephy, start.time = start.time, 
            end.time = end.time, node = node, nslices = nBins, 
            nodetype = nodetype)
    }
    if ("bamm-ratematrix" \%in\% class(ephy)) {
        if (!any(is.null(c(start.time, end.time, node)))) {
            stop("ERROR: You cannot specify start.time, end.time or node if the rate matrix is being provided. Please either provide the bammdata object instead or specify start.time, end.time or node in the creation of the bamm-ratematrix.")
        }
        rmat <- ephy
    }
    if (ratetype != "auto" & ratetype != "extinction" & ratetype != 
        "netdiv") {
        stop("ERROR: ratetype must be 'auto', 'extinction', or 'netdiv'.\n")
    }
    if (ephy$type == "trait" & ratetype != "auto") {
        stop("ERROR: If input object is of type 'trait', ratetype can only be 'auto'.")
    }
    if (ratetype == "auto" & ephy$type == "diversification") {
        rate <- rmat$lambda
        ratelabel <- "Speciation"
    }
    if (ratetype == "auto" & ephy$type == "trait") {
        rate <- rmat$beta
        ratelabel <- "trait rate"
    }
    if (ratetype == "extinction") {
        rate <- rmat$mu
        ratelabel <- "Extinction"
    }
    if (ratetype == "netdiv") {
        rate <- rmat$lambda - rmat$mu
        ratelabel <- "Net diversification"
    }
    maxTime <- max(rmat$times)
    if (!is.null(intervals)) {
        mm <- apply(rate, MARGIN = 2, quantile, intervals)
        poly <- list()
        q1 <- 1
        q2 <- nrow(mm)
        repeat {
            if (q1 >= q2) {
                break
            }
            a <- as.data.frame(cbind(rmat$times, mm[q1, ]))
            b <- as.data.frame(cbind(rmat$times, mm[q2, ]))
            b <- b[rev(rownames(b)), ]
            colnames(a) <- colnames(b) <- c("x", "y")
            poly[[q1]] <- rbind(a, b)
            q1 <- q1 + 1
            q2 <- q2 - 1
        }
    }
    if (!useMedian) {
        avg <- colMeans(rate)
    }
    else {
        avg <- unlist(apply(rate, 2, median))
    }
    if (smooth) {
        for (i in 1:length(poly)) {
            p <- poly[[i]]
            rows <- nrow(p)
            p[1:rows/2, 2] <- loess(p[1:rows/2, 2] ~ p[1:rows/2, 
                1], span = smoothParam)$fitted
            p[(rows/2):rows, 2] <- loess(p[(rows/2):rows, 2] ~ 
                p[(rows/2):rows, 1], span = smoothParam)$fitted
            poly[[i]] <- p
        }
        avg <- loess(avg ~ rmat$time, span = smoothParam)$fitted
    }
    if (plot) {
        if (!add) {
            plot.new()
            par(mar = mar)
            if (unique(xlim == "auto") & unique(ylim == "auto")) {
                plot.window(xlim = c(maxTime, 0), ylim = c(0, 
                  max(poly[[1]][, 2])))
                xMin <- maxTime
                xMax <- 0
                yMin <- 0
                yMax <- max(poly[[1]][, 2])
            }
            if (unique(xlim != "auto") & unique(ylim == "auto")) {
                plot.window(xlim = xlim, ylim = c(0, max(poly[[1]][, 
                  2])))
                xMin <- xlim[1]
                xMax <- xlim[2]
                yMin <- 0
                yMax <- max(poly[[1]][, 2])
            }
            if (unique(xlim == "auto") & unique(ylim != "auto")) {
                plot.window(xlim = c(maxTime, 0), ylim = ylim)
                xMin <- maxTime
                xMax <- 0
                yMin <- ylim[1]
                yMax <- ylim[2]
            }
            axis(at = c(1.3 * xMin, round(seq(xMin, xMax, length.out = xticks + 
                1))), labels = c(1.3 * xMin, round(seq(xMin, 
                xMax, length.out = xticks + 1))), cex.axis = cex.axis, 
                side = 1)
            axis(at = c(-0.2, seq(yMin, 1.2 * yMax, length.out = yticks + 
                1)), labels = c(-0.2, round(seq(yMin, 1.2 * yMax, 
                length.out = yticks + 1), digits = 1)), las = 1, 
                cex.axis = cex.axis, side = 2)
            mtext(side = 1, text = "Time since present", line = xline, 
                cex = cex)
            mtext(side = 2, text = ratelabel, line = yline, cex = cex)
        }
        if (!is.null(intervals)) {
            for (i in 1:length(poly)) {
                polygon(x = maxTime - poly[[i]][, 1], y = poly[[i]][, 
                  2], col = transparentColor(intervalCol, opacity), 
                  border = NA)
            }
        }
        lines(x = maxTime - rmat$time, y = avg, lwd = 3, col = avgCol)
    }
    else {
        return(list(poly = poly, avg = avg, times = rmat$time))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
