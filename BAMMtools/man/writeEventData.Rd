\name{writeEventData}
\alias{writeEventData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write a \code{bammdata} object to disk
}
\description{
\code{writEventData} takes a \code{bammdata} object and re-writes it back into a tree file and an event csv file.
}
\usage{
writeEventData(ephy, outtreefile, outeventfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ephy}{
\code{bammdata} object
}
  \item{outtreefile}{
the file name for outputting the tree
}  
  \item{outeventfile}{
the file name for outputting the event csv file
}
}
\author{
Huateng Huang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{subtreeBAMM}}
}
