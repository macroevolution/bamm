\name{generateControlFile}
\alias{generateControlFile}

\title{
	Generate control file for BAMM
}
\description{
	Generates a template diversification or trait control files for BAMM, while allowing the user to specify parameter values. 
}
\usage{
generateControlFile(file = "controlfile.txt", type = "diversification", params = NULL)
}
\arguments{
	\item{file}{
		destination file name with or without path
	}
	\item{type}{
		\code{diversification} or \code{trait}, depending on the desired BAMM analysis
	}
	\item{params}{
		list of parameters, see Details.
	}
}
\details{
	The user can supply parameters as a list, where the name of the list item is the name of the parameter as it appears in the control file, and the value of the list item is what will be placed in the contol file. 
	If a parameter is specified by the user, it will automatically be uncommented if it was commented in the template.
}

\references{
	\url{bamm-project.org}
}
\author{
	Pascal Title
}

\examples{
	
#Produce a blank template control file
generateControlFile(file = 'traitcontrol.txt', type='trait')

#Produce a customized control file
data(whales)

#get bamm priors to supply to control file
priors <- setBAMMpriors(whales, outfile = NULL)
rownames(priors) <- priors[,1]
priors[,2] <- round(priors[,2], 4)

generateControlFile(file = 'divcontrol.txt', params = list(
	treefile = 'whales.tre',
	globalSamplingFraction = '1',
	numberOfGenerations = '100000',
	overwrite = '1',
	lambdaInitPrior = priors['lambdaInitPrior', 2],
	lambdaShiftPrior = priors['lambdaShiftPrior', 2],
	muInitPrior = priors['muInitPrior', 2],
	poissonRatePrior = '1'))
}




