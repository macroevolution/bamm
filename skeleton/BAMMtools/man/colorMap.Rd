\name{colorMap}
\alias{colorMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
colorMap(x, pal, NCOLORS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, pal, NCOLORS) 
{
    dpal = c("BrBG", "PiYG", "PuOr", "RdBu", "RdGy", "RdYlBu", 
        "RdYlGn", "Spectral")
    colset = numeric(length(x))
    if (length(pal) == 3) {
        colpalette = colorRampPalette(pal, space = "Lab")(NCOLORS)
    }
    else if (pal \%in\% dpal) {
        colpalette = colorRampPalette(rev(brewer.pal(3, pal)), 
            space = "Lab")(NCOLORS)
    }
    else if (pal == "temperature") {
        colpalette = piratedColors(NCOLORS)
    }
    bks = quantile(x, seq(0, 1, length.out = (NCOLORS + 1)))
    for (i in 2:length(bks)) {
        if (i == 2) {
            colset[x < bks[2]] = colpalette[1]
        }
        else if (i == length(bks)) {
            colset[x >= bks[length(bks) - 1]] = colpalette[length(bks) - 
                1]
        }
        else {
            colset[x >= bks[i - 1] & x < bks[i]] = colpalette[i - 
                1]
        }
    }
    return(colset)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
