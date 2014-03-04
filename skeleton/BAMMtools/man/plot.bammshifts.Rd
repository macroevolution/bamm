\name{plot.bammshifts}
\alias{plot.bammshifts}

\title{
	Plot distinct rate shift configurations on a phylogeny
}
\description{
	Plots a random distinct rate shift configuration sampled by BAMM on a phylogeny.
}
\usage{
\method{plot}{bammshifts}(x, ephy, method="phylogram",pal="RdYlBu",
	 rank=NULL, index=NULL, spex="s", 
	 legend=TRUE, add.freq.text=TRUE,...)
}
\arguments{
	\item{x}{
		an object of class \code{bammshifts}.
	}
	\item{ephy}{
		an object of class \code{bammdata}.
	}
	\item{method}{
		a character string for which plotting method to use. 
		"phylogram" uses rectangular coordinates. "polar" uses
		polar coordinates.
	}
	\item{pal}{
		the color palette to use in \code{plot.bammdata}.
	}
	\item{rank}{
		the rank of the core shift configuration to plot. For the default (\code{NULL}) a random configuration is chosen. 
	}
	\item{index}{
		the posterior sample to plot. For the default (\code{NULL}) a random sample is chosen.
	}
	\item{spex}{
a character string indicating what type of macroevolutionary rates should be plotted. "s" (default) indicates
speciation rates, "e" indicates extinction rates, and any other character, e.g. "se", indicates net diversification
rates. Ignored if \code{ephy$type = "trait"}.
	}
	\item{legend}{
		a logical indicating whether to plot a legend.
	}
	\item{add.freq.text}{
		a logical indicating whether the frequency of each sampled shift configuration should be added to each plot.
	}
	\item{...}{
		other arguments to \code{plot.bammdata}.
	}	
}
\details{
A rate shift configuration is the set of nodes of the phylogeny where a shift occurs in the macroevolutionary rate dynamic of diversification or trait evolution. Each posterior sample is a potentially distinct rate shift configuration.
Different configurations may imply different macroevolutionary scenarios. This function helps visualize the different distinct rate shift configurations sampled by BAMM.

A core shift configuration is defined by a set of nodes that have a marginal shift probability greater than some threshold value that is specified in the function \code{distinctShiftConfigurations}. It is almost certain that more than one core shift configuration will be sampled by BAMM. Moreover, each core shift configuration may contain many subconfigurations. A subconfiguration contains the core shift node configuration and zero or more additional shift nodes that occur with low marginal probability.  

Points are added to the branches subtending the nodes of each rate configuration. The size of the point is proportional to the marginal probability that a shift occurs on a specific branch. If the instantaneous rate at a shift's origin represents an initial increase above the ancestral instantaneous rate the point is red. If the instantaneous rate at a shift's origin represents an initial decrease below the ancestral instantaneous rate the point is blue.   
	
}

\references{
	\url{bamm-project.org}
}
\author{
	Mike Grundler, Dan Rabosky
}
\seealso{
	\code{\link{distinctShiftConfigurations}}, \code{\link{plot.bammdata}}
}
\examples{
data(whales, events.whales, prior.whales)

ed <- getEventData(whales, events.whales, burnin=0.25)

# Here we get prior for identifying "significant" (core) rate shifts
priordist <- getBranchShiftPriors(whales, prior.whales)

sc <- distinctShiftConfigurations(ed, threshold=priordist)

plot(sc, ed)
}