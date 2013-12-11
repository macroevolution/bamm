


##########################################
#
#   computeJointMarginalShiftProbs
#
#   A tool for analyzing joint marginal distributions
#     of rate shift configurations from BAMM output.
#
#   args:
#	
#          ephy:              BAMM-data object
#          threshold:         WIll only consider 
#                              nodes with marginal shift probabilities
#                              greater than threshold
#                              Important to exclude singletons etc.
#  
#   Returns:
#          two matrices:
#          phi:               pairwise correlation matrix using phi coefficient
#          p.value            chi-square pvalue for the null hypothesis that there 
#                             is no significant association (+ or -) between the frequency
#                             of association between 2 shift nodes.
#                         
#
#
computeJointShiftCorrelations <- function(ephy, threshold=0.05){
	

	if ('bammdata' != class(ephy)){
		stop("Object ephy must be of class bammdata\n");
	}
	

	#First, restrict only to important nodes:
	#	those sampled at least with threshold frequency.
	
	rootnode <- length(ephy$tip.label) + 1;
	dd <- unlist(ephy$eventData);
	nodes <- dd[grep('node', names(dd))];
	nodes <- nodes[nodes != rootnode];
	tx <- table(nodes);
	tx <- tx / length(ephy$eventData);
	tx <- tx[tx >= threshold];	
	
	if (length(tx) == 0){
		stop("No nodes were sampled frequently enough to meet the specified threshold\n");
	}
	
	keepnodes <- as.numeric(names(tx));

	mm <- matrix(NA, nrow=length(keepnodes), ncol=length(keepnodes));

	# Make 4 matrices:
 
	n00 <- matrix(0, nrow=length(keepnodes), ncol=length(keepnodes));
	n11 <- n00;
	n01 <- n00;
	n10 <- n00;
	
	for (k in 1:length(ephy$eventData)){
		nodes <- ephy$eventData[[k]]$node;
		for (i in 1:(length(keepnodes)-1)){
			for (j in (i+1):length(keepnodes)){
				
				x <- keepnodes[i];
				y <- keepnodes[j];
				
				n00[i,j] <- n00[i,j] + (sum(c(x,y) %in% nodes) == 0);
				n11[i,j] <- n11[i,j] + (sum(c(x,y) %in% nodes) == 2);
				n10[i,j] <- n10[i,j] + ((x %in% nodes) & !(y %in% nodes));
				n01[i,j] <- n01[i,j] + ( !(x %in% nodes) & (y %in% nodes));						
			}
		}	
	}
	
	## Compute chi-square probability of each cell:
	
	chimat <- matrix(0, nrow=length(keepnodes), ncol=length(keepnodes));
	for (i in 1:(length(keepnodes)-1)){
		for (j in (i+1):length(keepnodes)){	
			tmp <- matrix(0, nrow=2, ncol=2);
			tmp[1,1] <- n00[i,j];
			tmp[1,2] <- n01[i,j];
			tmp[2,1] <- n10[i,j];
			tmp[2,2] <- n11[i,j];		
			chimat[i,j] <- chisq.test(tmp, simulate.p.value=T)$p.value;
		}	
	}
	
	phimat <- matrix(NA, nrow=length(keepnodes), ncol=length(keepnodes));
	
	for (i in 1:length(keepnodes)-1){
		for (j in (i+1):length(keepnodes)){
			
			x <- keepnodes[i];
			y <- keepnodes[j];
							
			num <- (n00[i,j]*n11[i,j] - n10[i,j]*n01[i,j]);
			denom <- (n00[i,j]+n01[i,j])*(n11[i,j]+n10[i,j])*(n00[i,j] + n10[i,j])*(n01[i,j] + n11[i,j]);
			phimat[i,j] <- num / sqrt(denom);
			phimat[j,i] <- phimat[i,j];
		}
	}
	diag(phimat) <- rep(1, nrow(phimat));
	rownames(phimat) <- keepnodes;
	colnames(phimat) <- keepnodes;

	
	chimat <- chimat + t(chimat);
	rownames(chimat) <- keepnodes;
	colnames(chimat) <- keepnodes;
	
	res <- list();
	res$phi <- phimat;
	res$p.value <- chimat;
	
	return(res);
	
}



