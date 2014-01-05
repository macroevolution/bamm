computeBayesFactors <-
function(postdata, priordata, burnin = 0.1, modelset = NULL, threshpost = 1, threshprior = 0, nbprior = FALSE, strict=FALSE){


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
		modelset <- unique(c(post[,2], prior[,2]));
	}else if (length(modelset) < 2){
		stop('\nInvalid modelset argument. This must be a vector of length > 1');
	}else{
		
	}	
	
	if (nbprior){
		
		subs <- prior[,2];
		if (length(subs) > 5000){
			subs <- subs[sample(1:length(subs), size=5000)];
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
