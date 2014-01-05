\name{computeBayesFactors}
\alias{computeBayesFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
temp
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeBayesFactors(postdata, priordata, burnin = 0.1, modelset = NULL, threshpost = 1, threshprior = 0, nbprior = FALSE, strict = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{postdata}{
%%     ~~Describe \code{postdata} here~~
}
  \item{priordata}{
%%     ~~Describe \code{priordata} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{modelset}{
%%     ~~Describe \code{modelset} here~~
}
  \item{threshpost}{
%%     ~~Describe \code{threshpost} here~~
}
  \item{threshprior}{
%%     ~~Describe \code{threshprior} here~~
}
  \item{nbprior}{
%%     ~~Describe \code{nbprior} here~~
}
  \item{strict}{
%%     ~~Describe \code{strict} here~~
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
function (postdata, priordata, burnin = 0.1, modelset = NULL, 
    threshpost = 1, threshprior = 0, nbprior = FALSE, strict = FALSE) 
{
    if (class(postdata) == "character") {
        post <- read.csv(postdata, header = T)
    }
    else if (class(postdata) == "data.frame") {
        post <- postdata
    }
    else {
        stop("invalid postdata argument (wrong class) in computeBayesFactors\n")
    }
    if (class(priordata) == "character") {
        prior <- read.csv(priordata, header = T)
    }
    else if (class(priordata) == "data.frame") {
        prior <- priordata
    }
    else {
        stop("invalid priordata argument (wrong class) in computeBayesFactors\n")
    }
    post <- post[floor(burnin * nrow(post)):nrow(post), ]
    prior <- prior[floor(burnin * nrow(prior)):nrow(prior), ]
    if (is.null(modelset)) {
        modelset <- unique(c(post[, 2], prior[, 2]))
    }
    else if (length(modelset) < 2) {
        stop("\nInvalid modelset argument. This must be a vector of length > 1")
    }
    else {
    }
    if (nbprior) {
        subs <- prior[, 2]
        if (length(subs) > 5000) {
            subs <- subs[sample(1:length(subs), size = 5000)]
        }
        resnb <- fitNegativeBinomial(subs)
        p1 <- resnb$size
        p2 <- resnb$mu
        tprior <- dnbinom(x = modelset, size = p1, mu = p2)
        names(tprior) <- as.character(modelset)
    }
    else {
        tprior <- table(prior[, 2])
    }
    tpost <- table(post[, 2])
    modelset <- sort(modelset)
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
            if (nbprior) {
                prior_odds <- priorf[i]/priorf[j]
            }
            else {
                prior_odds <- (priorf[i] + 1)/(priorf[j] + 1)
            }
            post_odds <- (postf[i] + 1)/(postf[j] + 1)
            ix <- modelset[i]
            ij <- modelset[j]
            if (!strict) {
                isGood1 <- (sum(post[, 2] == ix) >= threshpost) | 
                  (sum(post[, 2] == ij) >= threshpost)
                isGood2 <- (sum(prior[, 2] == ix) >= threshprior) | 
                  (sum(prior[, 2] == ij) >= threshprior)
                if (isGood1 & isGood2) {
                  mm[i, j] <- post_odds/prior_odds
                }
            }
            else {
                isGood1 <- (sum(post[, 2] == ix) >= threshpost) & 
                  (sum(post[, 2] == ij) >= threshpost)
                isGood2 <- (sum(prior[, 2] == ix) >= threshprior) & 
                  (sum(prior[, 2] == ij) >= threshprior)
                if (isGood1 & isGood2) {
                  mm[i, j] <- post_odds/prior_odds
                }
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
