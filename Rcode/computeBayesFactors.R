

#############################################################
#	computeBayesFactors
#
#
#   postfilename		=   MCMC output file from regular BAMM run
#				 			e.g., with sampleFromPriorOnly = 0
#			
#	priorfilename		=	MCMC output file from running BAMM with 
#                           sampleFromPriorOnly = 1
#
#
#	burnin				=	How many samples from posterior to discard
#							Uses fraction (e.g, 0.25 = 25% discarded)
#							Will also discard this same fraction from the prior.
#			
#	modelset			=	Integer set corresponding to models for which
#							you wish to compute pairwise Bayes Factors
#							e.g., 0:2 will compute all pairwise BF between models 
#							with 0 to 2 process 
#							(0 is a model with zero non-root processes)
#	
#							Will only compute Bayes Factors for the set of models
#							0:K that includes 99.5% of the sampled models. 
#
#   Returns:  matrix w pairwise Bayes Factors
#	By convention, the model with the higher index is the numerator for the calculation
#	e.g., M2 / M1 or M1 / M0, but never M0 / M1.
#
#	By default, odds ratios are computed as 
#			(prior_odds_M2  + 1   ) / (prior_odds_M1 + 1)
#		     where the 1 is added to both numerator and denominator 
#			 to avoid divide by zero erros 
#	
 
	
computeBayesFactors <- function(postfilename, priorfilename, burnin = 0.1, modelset = 0:5){

	if (length(modelset) < 2){
		stop('\nInvalid modelset argument. This must be a vector of length > 1');
	}
	
	
	post <- read.csv(postfilename, header=T);
	prior <- read.csv(priorfilename, header=T);
	
	post <- post[floor(burnin*nrow(post)):nrow(post), ];
	prior <- prior[floor(burnin*nrow(prior)):nrow(prior), ];

	
	tpost <- table(post$numevents);
	tprior <- table(prior$numevents);

	fprobs <- cumsum(tpost) / sum(tpost);
 	max_model <- NA;
 	if (length(fprobs) == 1){
 		max_model <- as.numeric(names(fprobs));
 	}else{
 		fprobs <- fprobs[fprobs < 0.995]; 		
 	 	max_model <- as.numeric(names(fprobs[length(fprobs)]));
	}
	
	if (max_model < max(modelset)){
		cat('*****************************************\n');
		cat('You have selected to compute Bayes Factors for models');
		cat('\n that were sampled very infrequently and for which\n');
		cat(' the Bayes Factors are likely to be (wildly ) inaccurate.\n');
		cat(' Consequently, the maximum rank of the models considered\n');
		cat(' will be constrained to <<< ', max_model, ' >>>\n');
		cat('*****************************************\n\n');

	}
	
	modelset <- modelset[modelset <= max_model];
	mset <- as.character(modelset);

	
	postf <- numeric(length(modelset));
	names(postf) <- mset;
	inboth <- intersect(mset, names(tpost));
	postf[inboth] <- tpost[inboth];
	
	priorf <- numeric(length(modelset));
	names(priorf) <- mset;
	inboth <- intersect(mset, names(tprior));
	priorf[inboth] <- tprior[inboth];
	
	mm <- matrix(NA, nrow=length(mset), ncol=length(mset));
	rownames(mm) <- mset;
	colnames(mm) <- mset;
	
	if (length(modelset) < 2){
		stop('\nError. Invalid model choice - is rank of specified model too high?\n');
	}
	
	for (i in 1:(length(modelset) - 1)){
		for (j in (i+1):length(modelset)){
					
			prior_odds <- (priorf[j] + 1) / (priorf[i] + 1);
			post_odds <- (postf[j] + 1) / (postf[i] + 1);
			
			mm[i , j] <- post_odds / prior_odds;
			
		}	
		
	}
	
	return(mm);
	
}


