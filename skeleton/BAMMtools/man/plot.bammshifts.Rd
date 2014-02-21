\name{plot.bammshifts}
\alias{plot.bammshifts}

\title{
	Plot distinct rate shift configurations on a phylogeny
}
\description{
	Plots each distinct rate shift configuration sampled by BAMM separately on a phylogeny.
}
\usage{
S3 method for class \code{bammshift}
plot(sc, ephy, plotmax=9, method="phylogram",pal="temperature",
	 spex="s",add.freq.text=TRUE, use.plot.bammdata=TRUE, 
	 send2pdf=FALSE)
}
\arguments{
	\item{sc}{
		an object of class \code{bammshifts}.
	}
	\item{ephy}{
		an object of class \code{bammdata}.
	}
	\item{plotmax}{
		an integer of the number of unique rate shift configurations to display.
	}
	\item{method}{
		a character string for which plotting method to use. 
		"phylogram" uses rectangular coordinates. "polar" uses
		polar coordinates.
	}
	\item{pal}{
		the color palette to use in \code{plot.bammdata}.
	}
	\item{spex}{
a character string indicating what type of macroevolutionary rates should be plotted. "s" (default) indicates
speciation rates, "e" indicates extinction rates, and any other character, e.g. "se", indicates net diversification
rates. Ignored if \code{ephy$type = "trait"}.
	}
	\item{add.freq.text}{
		a logical indicating whether the frequency of each sampled shift configuration should be added to each plot.
	}
	\item{use.plot.bammdata}{
		a logical indicating whether to use a phylorate phylogeny or not.
	}
	\item{send2pdf}{
		a logical indicating whether the output should be written to a PDF file. 
	}
}
\details{
A rate shift configuration is the set of nodes of the phylogeny where a shift occurs in the macroevolutionary rate dynamic of diversification or trait evolution. Each posterior sample is a potentially distinct rate shift configuration.
Different configurations may imply different macroevolutionary scenarios. This function helps visualize the different distinct rate shift configurations sampled by BAMM.

Points are added to the branches subtending the nodes of each rate configuration. The size of the point is proportional to the marginal probability that a shift occurs on a specific branch. If the instantaneous rate at a shift's origin represents an initial increase above the ancestral instantaneous rate the point is red. If the instantaneous rate at a shift's origin represents an initial decrease below the ancestral instantaneous rate the point is blue.   
	
}
\value{
	
}
\references{
	\link{bamm-project.org}
}
\author{
	Mike Grundler, Dan Rabosky
}
\seealso{
	\code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammdata}}
}
\examples{
library(BAMMtools)
data(whales, events.whales)
ed <- getEventData(whales, events.whales, burnin=0.25)

sc <- distinctShiftConfigurations(ed, threshold=0.05)
plot(sc, ed)
}