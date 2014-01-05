\name{piratedColors}
\alias{piratedColors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
piratedColors(n, palette = "temperature", alpha = 1, rgb = FALSE, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{palette}{
%%     ~~Describe \code{palette} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{rgb}{
%%     ~~Describe \code{rgb} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (n, palette = "temperature", alpha = 1, rgb = FALSE, 
    plot = FALSE) 
{
    if (n <= 0) 
        return(character(0))
    palette <- match.arg(palette, c("temperature", "blues"))
    x <- seq(0, 1, length = n)
    if (palette == "temperature") {
        r <- 1/(1 + exp(20 - 35 * x))
        g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
        b <- dnorm(x, 0.25, 0.15)/max(dnorm(x, 0.25, 0.15))
    }
    else {
        r <- 0.6 * x + 0.4 * x^2
        g <- 1.5 * x - 0.5 * x^2
        b <- 0.36 + 2.4 * x - 2 * x^2
        b[x > 0.4] <- 1
    }
    rgb.m <- matrix(c(r, g, b), ncol = 3, dimnames = list(NULL, 
        c("red", "green", "blue")))
    col <- mapply(rgb, r, g, b, alpha)
    if (rgb) 
        attr(col, "rgb") <- cbind(rgb.m, alpha)
    if (plot) {
        opar <- par("fig", "plt")
        par(fig = c(0, 1, 0, 0.7), plt = c(0.15, 0.9, 0.2, 0.95))
        plot(NA, xlim = c(-0.01, 1.01), ylim = c(-0.01, 1.01), 
            xlab = "Spectrum", ylab = "", xaxs = "i", yaxs = "i", 
            axes = FALSE)
        title(ylab = "Value", mgp = c(3.5, 0, 0))
        matlines(x, rgb.m, col = colnames(rgb.m), lty = 1, lwd = 3)
        matpoints(x, rgb.m, col = colnames(rgb.m), pch = 16)
        axis(1, at = 0:1)
        axis(2, at = 0:1, las = 1)
        par(fig = c(0, 1, 0.75, 0.9), plt = c(0.08, 0.97, 0, 
            1), new = TRUE)
        midpoints <- barplot(rep(1, n), col = col, border = FALSE, 
            space = FALSE, axes = FALSE)
        axis(1, at = midpoints, labels = 1:n, lty = 0, cex.axis = 0.6)
        par(opar)
    }
    return(col)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
