\name{samplingProbs}
\alias{samplingProbs}

\title{
Creates clade-specific sampling fractions
}
\description{
If the user would like to specify species sampling on a clade-by-clade basis, a sampling probability table can be provided to BAMM.
}
\usage{
samplingProbs(tree, cladeTable, cladeRichness = NULL, 
	globalSampling, output, writeToDisk = TRUE)
}

\arguments{
  \item{tree}{
an object of class \code{phylo}.
}
  \item{cladeTable}{
a dataframe with one column of species names and a second column of group assignment. 
}
  \item{cladeRichness}{
either \code{NULL} or a vector of species counts, named by clade names.
}
\item{globalSampling}{
percent sampling for the backbone of the phylogeny
}
\item{output}{
path + output file name
}
\item{writeToDisk}{
a logical, should the table be written to disk, defaults to \code{TRUE}.
}
}
\details{
This function handles two types of input: The cladeTable can either contain the species found in the phylogeny, 
along with clade assignment of those species, or it can contain more species than found in the phylogeny. 
If the table only contains those species in the phylogeny, then a vector \code{cladeRichness} must be provided 
that contains known clade richness. If the cladeTable contains more than the species in the phylogeny, 
then cladeRichness should be set to \code{NULL}.
The \code{globalSampling} value should represent the overall completeness of species sampling in terms of major clades.
See \url{bamm-project.org} for additional details.
}

\author{
Pascal Title
}
\value{
If \code{writeToDisk = TRUE}, then no object is returned. If \code{writeToDisk = FALSE}, then a dataframe is returned.
The resultant table must contain one row for each species in the phylogeny, along with clade assignment, and sampling fraction.
The first line must contain the overall sampling fraction for the phylogeny and must be written as tab-delimited, with no headers.
}


\examples{
#Generate dummy data
tree <- read.tree(text="(((t1:0.678555,(t2:0.064823,t3:0.064823):0.613732):0.291158,
	((t4:0.07812,t5:0.07812):0.019788,t6:0.097908):0.871805):0.030287,
	(t7:0.68478,(t8:0.259145,t9:0.259145):0.425635):0.31522);")
tree$tip.label <- paste(rep('Species',9),1:9,sep='')

spTable <- as.data.frame(matrix(nrow=9,ncol=2))
spTable[,1] <- tree$tip.label
spTable[1:3,2] <- 'cladeA'
spTable[4:6,2] <- 'cladeB'
spTable[7:9,2] <- 'cladeC'
richnessVec <- c(cladeA=5, cladeB=4, cladeC=12)

# Option 1: We have a table of clade assignment for the species in the tree, 
#	along with a vector of known clade richness
spTable
richnessVec
samplingProbs(tree, cladeTable = spTable, cladeRichness = richnessVec, 
	globalSampling = 1, writeToDisk = FALSE)

# Option 2: We have a table of known species, beyond the sampling in the phylogeny
spTable <- rbind(spTable, c('Species10','cladeA'),c('Species11','cladeA'),
	c('Species12','cladeC'), c('Species13','cladeC'),c('Species14','cladeC'))

spTable

samplingProbs(tree, cladeTable = spTable, cladeRichness = NULL, 
	globalSampling = 0.9, writeToDisk = FALSE)

}

