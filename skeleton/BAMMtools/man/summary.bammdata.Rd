\name{summary.bammdata}
\alias{summary.bammdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp title
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.bammdata(ephy, opt = "product", display = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
%%     ~~Describe \code{ephy} here~~
}
  \item{opt}{
%%     ~~Describe \code{opt} here~~
}
  \item{display}{
%%     ~~Describe \code{display} here~~
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
function (ephy, opt = "product", display = 10) 
{
    cat("\nAnalyzed", length(ephy$eventData), "posterior samples\n")
    shiftsindex = maximumShiftCredibilityTree(ephy, maximize = opt)$sampleindex
    shiftnodes = getShiftNodesFromIndex(ephy, shiftsindex)
    if (length(shiftnodes) > 1) {
        cat("\nMaximum shift credibility tree has", length(shiftnodes), 
            "shifts\n\n")
        cat("Shifts occur leading to nodes:", shiftnodes, "\n\n")
        cat("Optimality type:", opt, "\n\n")
    }
    else if (length(shiftnodes) == 1) {
        cat("\nMaximum shift credibility tree has 1 shift\n\n")
        cat("Shift occurs leading to node:", shiftnodes, "\n\n")
        cat("Optimality type:", opt, "\n\n")
    }
    else {
        cat("\nMaximum shift credibility tree has 0 shifts\n\n")
        cat("Optimality type:", opt, "\n\n")
    }
    cat("Shift posterior distribution:\n\n")
    fev = sapply(ephy$eventData, nrow)
    disp = tabulate(fev)
    disp = disp/sum(disp)
    disp = data.frame(cbind(seq.int(0, length(disp) - 1, 1), 
        signif(disp, 2)))
    if (nrow(disp) <= display) {
        write.table(format(disp, justify = "left", width = 10), 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
    }
    else {
        if (nrow(disp) - display == 1) {
            write.table(format(disp[1:display, ], justify = "left", 
                width = 10), col.names = FALSE, row.names = FALSE, 
                quote = FALSE)
            cat("... omitted 1 row\n\n")
        }
        else {
            write.table(format(disp[1:display, ], justify = "left", 
                width = 10), col.names = FALSE, row.names = FALSE, 
                quote = FALSE)
            cat("... omitted", nrow(disp) - display, "rows\n\n")
        }
    }
    invisible(list(posterior = tabulate(fev), mscshiftnodes = shiftnodes))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
