#############################################################
#	computeBayesFactors
#
#
#   postdata		=   MCMC output file from regular BAMM run
#				 			e.g., with sampleFromPriorOnly = 0
#							OR a dataframe			
#
#	priordata		=	MCMC output file from running BAMM with 
#                           sampleFromPriorOnly = 1
#							OR a dataframe
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
#							If is.null(modelset), this will assume modelset consists of 
#							all sampled models
#	
#
#  threshpost, threshprior	=   Will only compute BF for a model comparison where 
#	 						at least one of the models has been sampled at least
#							thresh times. This avoids comparisons between two models
#							that were very rarely or never sampled, which always implies
#							highly inaccurate posterior or prior probabilities	
#	nbprior				=   use negative binomial distribution to 
#								approximate the full prior distribution   
#							This runs into trouble in some cases. When the prior is approximated
#								with a high level of accuracy, but the posterior odds are poorly estimated
#								this tends to fail.
#	strict 		        =   logical. If TRUE, requires that both 
#							models i and j be sampled at least threshpost or threshprior times.
#
#
#   Returns:  matrix w pairwise Bayes Factors
#			  BF_{i, j} is the Bayes factor between model i (numerator)
#							and model j (denominator)
#
#	By default, odds ratios are computed as 
#			(prior_odds_M2  + 1   ) / (prior_odds_M1 + 1)
#		     where the 1 is added to both numerator and denominator 
#			 to avoid divide by zero erros 
#	
#   Dependency on BAMM MCMC output: if order of output columns 
#		changes, it will break this function.
# 
#   This function can be very difficult to use.

computeBayesFactors <- function(postdata, priordata, burnin = 0.1, modelset = NULL, threshpost = 1, threshprior = 0, nbprior = FALSE, strict=FALSE){


	if (class(postdata) == 'character'){
		post <- read.csv(postdata, header=T);
	}else if (class(postdata) == 'data.frame'){
		post <- postdata;
	}else{
		stop("invalid postdata argument (wrong class) in computeBayesFactors\n");
	}

	if (class(priordata) == 'character'){
		prior <- read.csv(priordata, header=T);
	}else if (class(priordata) == 'data.frame'){
		prior <- priordata;
	}else{
		stop("invalid priordata argument (wrong class) in computeBayesFactors\n");
	}
 
	post <- post[floor(burnin*nrow(post)):nrow(post), ];
	prior <- prior[floor(burnin*nrow(prior)):nrow(prior), ];
 
	
	if (is.null(modelset)){
		modelset <- unique(post[,2]);
	}else if (length(modelset) < 2){
		stop('\nInvalid modelset argument. This must be a vector of length > 1');
	}else{
		
	}	
	
	if (nbprior){
		
		subs <- prior[,2];
		if (length(subs) > 5000){
			subs <- subs(sample(1:length(subs), size=5000));
		}

		resnb <- fitNegativeBinomial(subs);

		p1 <- resnb$size;
		p2 <- resnb$mu;

		tprior <- dnbinom(x=modelset, size=p1, mu=p2);
		names(tprior) <- as.character(modelset);
		
	}else{
		tprior <- table(prior[,2]);	
	}



	tpost <- table(post[,2]);
	
	
	modelset <- sort(modelset);
	
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
	
 
	
	
	for (i in 1:length(modelset)){
		
		for (j in 1:length(modelset)){
					
			if (nbprior){
				prior_odds <- priorf[i] / priorf[j];				
			}else{
				prior_odds <- (priorf[i] + 1) / (priorf[j] + 1);				
			}
			post_odds <- (postf[i] + 1) / (postf[j] + 1);
			
			ix <- modelset[i];
			ij <- modelset[j];
 
			
			if (!strict){
				isGood1 <- (sum(post[,2] == ix) >= threshpost) | (sum(post[,2] == ij) >= threshpost); 
				isGood2 <- (sum(prior[,2] == ix) >= threshprior) | (sum(prior[,2] == ij) >= threshprior); 
		
				if (isGood1 & isGood2){
					mm[i , j] <- post_odds / prior_odds;	
				}				
			}else{
				# All models must be sampled at least once.
				isGood1 <- (sum(post[,2] == ix) >= threshpost) & (sum(post[,2] == ij) >= threshpost); 
				isGood2 <- (sum(prior[,2] == ix) >= threshprior) & (sum(prior[,2] == ij) >= threshprior); 				
				
				if (isGood1 & isGood2){
					mm[i , j] <- post_odds / prior_odds;	
				}					
			}
			

			
		}	
		
	}

	
	return(mm);
	
}



