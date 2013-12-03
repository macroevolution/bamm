\name{computeBayesFactors}
\alias{computeBayesFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeBayesFactors(postfilename, priorfilename, burnin = 0.1, modelset = NULL, threshold = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{postfilename}{
%%     ~~Describe \code{postfilename} here~~
}
  \item{priorfilename}{
%%     ~~Describe \code{priorfilename} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{modelset}{
%%     ~~Describe \code{modelset} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
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
function (postfilename, priorfilename, burnin = 0.1, modelset = NULL, 
    threshold = 1) 
{
    post <- read.csv(postfilename, header = T)
    prior <- read.csv(priorfilename, header = T)
    post <- post[floor(burnin * nrow(post)):nrow(post), ]
    prior <- prior[floor(burnin * nrow(prior)):nrow(prior), ]
    tpost <- table(post[, 2])
    tprior <- table(prior[, 2])
    if (is.null(modelset)) {
        modelset <- unique(post[, 2])
    }
    else if (length(modelset) < 2) {
        stop("\nInvalid modelset argument. This must be a vector of length > 1")
    }
    else {
    }
    mset <- as.character(modelset)
    postf <- numeric(length(modelset))
    names(postf) <- mset
    inboth <- intersect(mset, names(tpost))
    postf[inboth] <- tpost[inboth]
    priorf <- numeric(length(modelset))
    names(priorf) <- mset
    inboth <- intersect(mset, names(tprior))
    priorf[inboth] <- tprior[inboth]
    mm <- matrix(NA, nrow = length(mset), ncol = length(mset))
    rownames(mm) <- mset
    colnames(mm) <- mset
    if (length(modelset) < 2) {
        stop("\nError. Invalid model choice - is rank of specified model too high?\n")
    }
    for (i in 1:length(modelset)) {
        for (j in 1:length(modelset)) {
            prior_odds <- (priorf[i] + 1)/(priorf[j] + 1)
            post_odds <- (postf[i] + 1)/(postf[j] + 1)
            ix <- modelset[i]
            ij <- modelset[j]
            isGood_i <- sum(post[, 2] == ix) >= threshold
            isGood_j <- sum(post[, 2] == ij) >= threshold
            if (isGood_i | isGood_j) {
                mm[i, j] <- post_odds/prior_odds
            }
        }
    }
    return(mm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
