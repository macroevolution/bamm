computeBayesFactors <-
function(postfilename, priorfilename, burnin = 0.1, modelset = NULL, threshold = 1){


	
	post <- read.csv(postfilename, header=T);
	prior <- read.csv(priorfilename, header=T);
	
	post <- post[floor(burnin*nrow(post)):nrow(post), ];
	prior <- prior[floor(burnin*nrow(prior)):nrow(prior), ];

	tpost <- table(post[,2]);
	tprior <- table(prior[,2]);
	
	
	if (is.null(modelset)){
		modelset <- unique(post[,2]);
	}else if (length(modelset) < 2){
		stop('\nInvalid modelset argument. This must be a vector of length > 1');
	}else{
		
	}
	
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
					
			prior_odds <- (priorf[i] + 1) / (priorf[j] + 1);
			post_odds <- (postf[i] + 1) / (postf[j] + 1);
			
			ix <- modelset[i];
			ij <- modelset[j];
			isGood_i <- sum(post[,2] == ix) >= threshold;
			isGood_j <- sum(post[,2] == ij) >= threshold;
		
			if (isGood_i | isGood_j){
				mm[i , j] <- post_odds / prior_odds;	
			}
			
		}	
		
	}

	
	# for (i in 1:(length(modelset) - 1)){
		# for (j in (i+1):length(modelset)){
					
			# prior_odds <- (priorf[j] + 1) / (priorf[i] + 1);
			# post_odds <- (postf[j] + 1) / (postf[i] + 1);
			
			# mm[i , j] <- post_odds / prior_odds;
			
		# }	
		
	# }
	
	return(mm);
	
}
