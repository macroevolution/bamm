#############################################################
#
#	getShiftNodesFromIndex (....)
#
#	Args: ephy	=	object of class 'bamm-data'
#		  index =   the index of the sample you wish to view, e.g.,
#					 if index = 5, this will give you the nodes subtending
#                    all branches with rate shifts for the 5th sample
#					 from the posterior in your BAMM-data object.								 
#
#	
#	Returns: 		- a vector of the nodes where shifts occurred, excluding the root.
#					Note: if NO shifts occured, this will return a 
#							numeric object of length zero
# 

getShiftNodesFromIndex <- function(ephy, index){
	
	if (index > length(ephy$eventData)){
		cat("Error. Attempting to access non-existent element from BAMM-data object\n");
		cat("You have << ", length(ephy$eventData), " >>> samples in your BAMM-data object\n");
		cat("Value of index must be no greater than the number of samples\n\n");
		stop();
	}
	root <- length(ephy$tip.label) + 1;
	nodes <- ephy$eventData[[index]]$node;
	nodes <- nodes[nodes != root];

	return(nodes);
}


as.phylo.bammdata <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}		
	
	newphylo <- list();
	newphylo$edge <- ephy$edge;
	newphylo$Nnode <- ephy$Nnode;
	newphylo$tip.label <- ephy$tip.label;
	newphylo$edge.length <- ephy$edge.length;
	class(newphylo) <- 'phylo';
	return(newphylo);
}


#	type		= the type of analysis to autotune (speciationextinction or trait)
#					should probably correspond to default bamm modeltype options
#	burn		=	fraction of samples from the tuning output file to discard 
#						(recommend 0.1)
#	outfile		= 	file where summary will be stored.
#	
#	opt_prob	= theoretical optimum mcmc acceptance rate
#				0.234 is probably too low, but 
#				anything between 0.23 and 0.4 is probably fine.
#				(0.234 being optimum for infinite-dimensional target space)
#	
#
#
# Parameters by index, for diversification:
#		0 = change number of events (not tuned)
#		1 = event position scale (event location)
#		2 = updateEventRateScale
#		3 = updateLambdaInitScale
#		4 = updateLambdaShiftScale
#		5 = updateMuInitScale
#	

autotune_BAMM_MCMC <- function(fname, type = 'speciationextinction', burn=0.1, opt_prob = 0.234, outfile = 'Autotune_results.txt', min.set = 0.01){
	
	xx <- read.csv(fname, header=F);
	mn <- floor(0.1*nrow(xx));
	xx <- xx[mn:nrow(xx), ];
	
 	if (type == 'speciationextinction'){
 		
 		parnames = c('updateEventLocationScale', 'updateEventRateScale', 
 				'updateLambdaInitScale', 'updateLambdaShiftScale', 'updateMuInitScale');
 		
 		opt_vals <- numeric(length(parnames));
 		comments <- character(length(parnames));
 		
 		for (i in 1:5){
 			xtemp <- xx[xx[,1] == i, ];
 			fit_temp <- glm(xtemp[,3] ~ xtemp[,2], family = 'binomial');
 			opt_vals[i] <- (logit(opt_prob) - fit_temp$coefficients[1]) / fit_temp$coefficients[2];
 			if (opt_vals[i] < 0){
 				comments[i] <- ' Warning: this parameter not tuned. Setting to minimum value';
 				opt_vals[i] <- min.set;
 			}
 		
 		}		
 		
 		setup_outfile(outfile, parnames, opt_vals, comments);	
 		
 	}else if (type == 'trait'){
 		
 		cat("trait auto-tune not yet supported\n");
 	}
	
	
	
}


setup_outfile <- function(fname, parnames, parvalues, comments){
	cat('#\n#\tFile generated ', as.character(Sys.time()), '\n', file = fname, append=F);
	cat('#\t\t\tby BAMMtools::autotune_BAMM_MCMC\n', file = fname, append=T);	
	cat('#\n#\tThis file contains autotuned parameters\n', file = fname, append=T);
	cat('#\tfor a BAMM analysis\n#\n', file = fname, append=T);
	cat('#\tCopy the block of parameter below\n', file = fname, append=T);	
	cat('#\tinto the MCMC_OPERATORS section of your controlfile\n', file = fname, append=T);		
	cat('#\tEnsure that the default MCMC_SCALING_OPERATORS section is deleted\n', file = fname, append=T);	
	cat('#\tstart MCMC_SCALING_OPERATORS block <<<<\n\n', file = fname, append=T);	
	for (i in 1:length(parnames)){
		cat(parnames[i], ' = ', parvalues[i], sep='', file = fname, append=T);
		if (comments[i] != ''){
			cat('   #', comments[i], sep='', file = fname, append=T);				
		}
		cat('\n', file = fname, append=T);		
	}
	cat('\n#\t>>>>>end MCMC_SCALING_OPERATORS block', file = fname, append=T);	
	cat('\n#\tNote: if you received a warning message after a parameter,', file = fname, append=T);	
	cat('\n#\tthis just means that the auto-tuning could not identify a best target value', file = fname, append=T);	
	cat('\n#\tfor this scaling operator. You do not need to worry about this.', file = fname, append=T);	
	cat('\n#\tHowever, if MCMC is performing poorly, you may consider manuually', file = fname, append=T);	
	cat('\n#\tIncreasing or decreasing the value of the parameter', file = fname, append=T);	
}

inv.logit <- function (x, min = 0, max = 1) 
{
    p <- exp(x)/(1 + exp(x))
    p <- ifelse(is.na(p) & !is.na(x), 1, p)
    p * (max - min) + min
}

logit <- function (x, min = 0, max = 1) 
{
    p <- (x - min)/(max - min)
    log(p/(1 - p))
}
