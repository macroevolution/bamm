\name{getEventData}
\alias{getEventData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read event data from BAMM output
}
\description{
This function reads in BAMM event data and converts it into a bamm-data object. Most other functions in BammTools require this data object.
}
\usage{
getEventData(phy, eventfilename, burnin = 0, nsamples = NULL, verbose = FALSE, assign.type = "new_way", type = "diversification", header = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
The model tree (time calibrated tree analyzed with BAMM).
}
  \item{eventfilename}{
file where event data are stored, in BAMM output format.
}
  \item{burnin}{
number of samples from posterior to discard, as fraction (e.g, 0.25 = 25\% discarded).
}
  \item{nsamples}{
number of samples from posterior to include.If NULL, includes ALL samples (could take long time to run).
}
  \item{verbose}{
verbose print output for bug tracking.
}
  \item{assign.type}{
ignore for now.
}
  \item{type}{
"diversification" if output is from speciationextinction BAMM, "traits" if output is from BAMM traits. 
}
  \item{header}{
boolean to flag whether eventfilename contains a header.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{edge}{Edge matrix of model tree}
\item{Nnode}{Number of nodes in model tree}
\item{tip.label}{}

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
function (phy, eventfilename, burnin = 0, nsamples = NULL, verbose = FALSE, 
    assign.type = "new_way", type = "diversification", header = TRUE) 
{
    if (type != "diversification" & type != "traits") {
        stop("Invalid 'type' specification. Should be 'diversification' or 'traits'")
    }
    if (any(is.null(c(phy$begin, phy$end)))) {
        phy <- getStartStopTimes(phy)
    }
    bt <- branching.times(phy)
    eventVectors <- list()
    eventData <- list()
    tipStates <- list()
    eventBranchSegs <- list()
    tipLambda <- list()
    cat("Reading event datafile: ", eventfilename, "\n\t\t...........")
    x <- read.csv(eventfilename, header = header, stringsAsFactors = FALSE)
    uniquegens <- sort(unique(x[, 1]))
    cat("\nRead a total of ", length(uniquegens), " samples from posterior\n")
    samplestart <- uniquegens[floor(burnin * length(uniquegens))]
    if (!length(samplestart)) {
        samplestart <- 0
    }
    uniquegens <- uniquegens[uniquegens >= samplestart]
    if (is.null(nsamples)) {
        nsamples <- length(uniquegens)
    }
    else if (nsamples > length(uniquegens)) {
        nsamples <- length(uniquegens)
    }
    goodsamples <- uniquegens[seq.int(1, length(uniquegens), 
        length.out = nsamples)]
    cat("\nDiscarded as burnin: GENERATIONS < ", goodsamples[1])
    cat("\nAnalyzing ", length(goodsamples), " samples from posterior\n")
    numberEvents <- length(goodsamples)
    cat("\nSetting recursive sequence on tree...\n")
    if (assign.type == "new_way") 
        phy <- getRecursiveSequence(phy)
    cat("\nDone with recursive sequence\n\n")
    meanTipMu <- numeric(length(phy$tip.label))
    meanTipLambda <- numeric(length(phy$tip.label))
    for (i in 1:length(goodsamples)) {
        tmpEvents <- x[x[, 1] == uniquegens[i], ]
        if (verbose) 
            cat("Processing event: ", i, "\n")
        t1 <- tmpEvents[, 2]
        t2 <- tmpEvents[, 3]
        tm <- tmpEvents[, 4]
        lam1 <- tmpEvents[, 5]
        lam2 <- tmpEvents[, 6]
        if (type == "diversification") {
            mu1 <- tmpEvents[, 7]
            mu2 <- tmpEvents[, 8]
        }
        else {
            mu1 <- rep(0, nrow(tmpEvents))
            mu2 <- rep(0, nrow(tmpEvents))
        }
        tipMu <- list()
        nodeVec <- numeric(nrow(tmpEvents))
        for (k in 1:length(nodeVec)) {
            if (is.na(t2[k])) {
                nodeVec[k] <- which(phy$tip.label == t1[k])
            }
            else {
                tipnodes <- c(which(phy$tip.label == t1[k]), 
                  which(phy$tip.label == t2[k]))
                nodeVec[k] <- getMRCA(phy, tipnodes)
            }
        }
        if (sum(nodeVec == 0) > 0) {
            stop("Failed to assign event to node\n")
        }
        dftemp <- data.frame(node = nodeVec, time = tm, lam1 = lam1, 
            lam2 = lam2, mu1 = mu1, mu2 = mu2, stringsAsFactors = FALSE)
        dftemp <- dftemp[order(dftemp$time), ]
        dftemp$index <- 1:nrow(dftemp)
        rownames(dftemp) <- NULL
        tphy <- phy
        tphy$statevec <- rep(1, nrow(tphy$edge))
        if (assign.type == "old_way") {
            for (k in 1:nrow(dftemp)) {
                tphy <- recursivelySetNodeStates(tphy, node = dftemp$node[k], 
                  state = k)
            }
        }
        else if (assign.type == "new_way") {
            if (nrow(dftemp) > 1) {
                for (k in 2:nrow(dftemp)) {
                  s1 <- which(tphy$downseq == dftemp$node[k])
                  s2 <- which(tphy$downseq == tphy$lastvisit[dftemp$node[k]])
                  descSet <- tphy$downseq[s1:s2]
                  isDescendantNode <- tphy$edge[, 2] \%in\% descSet
                  tphy$statevec[isDescendantNode] <- k
                }
            }
        }
        tmpEventSegMat <- matrix(0, nrow = (max(phy$edge) + nrow(dftemp) - 
            2), ncol = 4)
        non.root <- c(1:length(phy$tip.label), (length(phy$tip.label) + 
            2):max(phy$edge))
        pos <- 1
        for (k in non.root) {
            events.on.branch <- dftemp[dftemp$node == k, ]
            if (nrow(events.on.branch) == 0) {
                tmpEventSegMat[pos, 1] <- k
                tmpEventSegMat[pos, 2] <- tphy$begin[tphy$edge[, 
                  2] == k]
                tmpEventSegMat[pos, 3] <- tphy$end[tphy$edge[, 
                  2] == k]
                tmpEventSegMat[pos, 4] <- tphy$statevec[tphy$edge[, 
                  2] == k]
            }
            else {
                events.on.branch <- events.on.branch[order(events.on.branch$time), 
                  ]
                fBranch <- phy$edge[, 2] == k
                start.times <- c(phy$begin[fBranch], events.on.branch$time)
                stop.times <- c(events.on.branch$time, phy$end[fBranch])
                parent <- phy$edge[, 1][phy$edge[, 2] == k]
                if (parent == (length(phy$tip.label) + 1)) {
                  proc.set <- c(1, events.on.branch$index)
                }
                else {
                  proc.set <- c(tphy$statevec[tphy$edge[, 2] == 
                    parent], events.on.branch$index)
                }
                zzindex <- pos:(pos + nrow(events.on.branch))
                tmpEventSegMat[zzindex, 1] <- rep(k, length(zzindex))
                tmpEventSegMat[zzindex, 2] <- start.times
                tmpEventSegMat[zzindex, 3] <- stop.times
                tmpEventSegMat[zzindex, 4] <- proc.set
            }
            pos <- pos + 1 + nrow(events.on.branch)
        }
        eventBranchSegs[[i]] <- tmpEventSegMat
        tipstates <- numeric(length(phy$tip.label))
        for (k in 1:length(tipstates)) {
            tipstates[k] <- tphy$statevec[tphy$edge[, 2] == k]
        }
        stoptime <- max(branching.times(phy))
        tiplam <- dftemp$lam1[tipstates] * exp(dftemp$lam2[tipstates] * 
            (stoptime - dftemp$time[tipstates]))
        tipmu <- dftemp$mu1[tipstates]
        meanTipMu <- meanTipMu + tipmu/nsamples
        meanTipLambda <- meanTipLambda + tiplam/nsamples
        eventData[[i]] <- dftemp
        eventVectors[[i]] <- tphy$statevec
        numberEvents[i] <- nrow(dftemp)
        tipStates[[i]] <- tipstates
        tipLambda[[i]] <- tiplam
        tipMu[[i]] <- tipmu
    }
    phy$numberEvents <- numberEvents
    phy$eventData <- eventData
    phy$eventVectors <- eventVectors
    phy$tipStates <- tipStates
    phy$tipLambda <- tipLambda
    phy$meanTipLambda <- meanTipLambda
    phy$eventBranchSegs <- eventBranchSegs
    phy$tipMu <- tipMu
    phy$meanTipMu <- meanTipMu
    if (type == "diversification") {
        phy$type = "diversification"
    }
    else {
        phy$type = "traits"
    }
    class(phy) <- c("phylo", "bamm-data")
    return(phy)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
