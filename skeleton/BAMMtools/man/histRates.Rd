\name{histRates}
\alias{histRates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
histRates(rates, pal, NCOLORS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rates}{
%%     ~~Describe \code{rates} here~~
}
  \item{pal}{
%%     ~~Describe \code{pal} here~~
}
  \item{NCOLORS}{
%%     ~~Describe \code{NCOLORS} here~~
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
function (rates, pal, NCOLORS) 
{
    opar = par(no.readonly = TRUE)
    fx = density(rates)
    dpal = c("BrBG", "PiYG", "PuOr", "RdBu", "RdGy", "RdYlBu", 
        "RdYlGn", "Spectral")
    if (length(pal) == 3) {
        rate.colors = colorRampPalette(pal, space = "Lab")(NCOLORS)
    }
    else if (pal \%in\% dpal) {
        rate.colors = colorRampPalette(rev(brewer.pal(3, pal)), 
            space = "Lab")(NCOLORS)
    }
    else if (pal == "temperature") {
        rate.colors = piratedColors(NCOLORS)
    }
    qx = quantile(rates, seq(0, 1, length.out = NCOLORS + 1))
    cat("Click once where you want the lower left corner of the figure\n")
    cxy = locator(n = 1)
    xc = (grconvertX(cxy$x, to = "ndc"))
    yc = (grconvertY(cxy$y, to = "ndc"))
    ofs = min(1 - xc, 1 - yc)
    fig = c(xc, xc + ofs, yc, yc + ofs)
    par(fig = fig, new = TRUE, xpd = TRUE, mar = c(1.5, 1.5, 
        0.25, 0.25))
    plot.new()
    plot.window(xlim = c(0, max(fx$x)), ylim = c(0, max(fx$y)))
    for (i in 1:length(fx$x)) {
        index = which(qx > fx$x[i])[1]
        if (is.na(index)) 
            break
        if (index > 1) 
            index = index - 1
        bcol = rate.colors[index]
        segments(fx$x[i], fx$y[i], fx$x[i], 0, lend = 2, col = bcol)
    }
    axis(1, signif(seq(0, max(fx$x), length.out = 5), 2), pos = 0, 
        cex.axis = 0.75, tcl = NA, mgp = c(0, 0.25, 0))
    axis(2, round(seq(0, max(fx$y), length.out = 3), 0), las = 1, 
        pos = 0, cex.axis = 0.75, tcl = NA, mgp = c(0, 0.25, 
            0))
    mtext("Evolutionary Rate", 1, line = 1, cex = 0.75)
    mtext("Density", 2, line = 1, cex = 0.75)
    par(opar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
