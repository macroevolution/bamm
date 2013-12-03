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
