#############################################################
#
#	getMeanBranchLengthTree(....)
#
#	Should be capable of computing this from 2 different types of data:
#		1. The bamm-data object
#		2. A list of phylogenetic trees,
#			e.g., as you would get directly as output from BAMM with the 
#			marginal branch length trees.
#		Arg: (for now) 
#			obj	= a list of multiple phylogenetic trees, where each tree
#						has branch length equal to the mean rate on that branch
#						OR a bamm-data object
#			ndr = TRUE or FALSE. if TRUE the mean branch length tree will have the
#				mean net rate on each branch, i.e., lambda - mu. if FALSE the mean
#				rate will simply be lambda.  For trait data ndr will make no difference  

getMeanBranchLengthTree <- function(obj,ndr=TRUE){
	
	if('bamm-data' %in% class(obj))
	{
		#v <- list(edge = obj$edge, Nnode = obj$Nnode,tip.label = obj$tip.label, edge.length = obj$edge.length);
		#attributes(v) <- list(names = c("edge","Nnode","tip.label","edge.length"),class="phylo",order="cladewise");
		v <- as.phylo.bammdata(obj);
		
		if(ndr)
		{
			if(obj$type == 'diversification'){
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$lambda_branch_matrix) - rowMeans(obj$mu_branch_matrix);
			}
			else if(obj$type == 'traits'){
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$beta_branch_matrix);
			}
		}
		else
		{
			if(obj$type == 'diversification'){
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$lambda_branch_matrix);
			}
			else if(obj$type == 'traits'){
				obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
				el  <- rowMeans(obj$beta_branch_matrix);
			}	
		}
		v$edge.length <- el;
		obj <- list();
		obj$phy <- v;
		obj$median <- median(el);
		obj$mean <- mean(el);			
	}
	else if(class(obj[[1]]) == 'phylo')
	{
		el = rowMeans(sapply(obj,with,edge.length));

		#v <- obj[[1]];
		#meanvec <- numeric(length(obj));
 
		#el <- v$edge.length;
		#meanvec[1] <- mean(v$edge.length);
		#for (i in 2:length(obj)){
		#	el <- el + obj[[i]]$edge.length;
		#	meanvec[i] <- mean(obj[[i]]$edge.length);
		#}
		#el <- el / length(obj);
		#v$edge.length <- el;
		
		obj <- list();
		obj$phy <- obj[[1]];
		obj$median <- median(el);
		obj$mean <- mean(el);
	}
	else
	{
		stop("Method not implemented for supplied object class");
	}
	return(obj);		
}
