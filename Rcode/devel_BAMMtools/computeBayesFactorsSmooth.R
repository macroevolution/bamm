
# Efficient Bayes factor estimation 

computeBayesFactorsSmooth <- function(postdata, priordata, burnin = 0.1, modelset = NULL, threshold = 1, bw=5){


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
	
	modelset <- sort(modelset);
	post_kd <- density(post[,2], from=min(modelset), to=max(modelset), bw=bw);
	prior_kd <- density(prior[,2], from=min(modelset), to=max(modelset), bw=bw);
	
	mm <- matrix(NA, nrow=length(modelset), ncol=length(modelset));
	rownames(mm) <- modelset;
	colnames(mm) <- modelset;
	
	if (length(modelset) < 2){
		stop('\nError. Invalid model choice - is rank of specified model too high?\n');
	}
	
 	indices <- numeric(length(modelset));
 	for (i in 1:length(modelset)){
 		kk <- abs(post_kd$x - modelset[i]);
 		indices[i] <- which(kk == min(kk))[1];
 	}
	
	
	for (i in 1:length(modelset)){
		
		for (j in 1:length(modelset)){
			
		 
			prior_odds <- (prior_kd$y[indices[i]] / prior_kd$y[indices[j]]);
			post_odds <- post_kd$y[indices[i]] / post_kd$y[indices[j]];

			ix <- modelset[i];
			ij <- modelset[j];
			isGood_i <- sum(post[,2] == ix) >= threshold;
			isGood_j <- sum(post[,2] == ij) >= threshold;
			
 			if (isGood_i | isGood_j){
				mm[i , j] <- post_odds / prior_odds;	
			}
			
		}	
		
	}

	
	return(mm);
	
}