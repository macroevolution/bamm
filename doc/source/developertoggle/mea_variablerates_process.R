
# Functions to create matrices 
# of branch and tip-specific rates from Moore et al (2016) input files
# as available on Dryad, http://datadryad.org/resource/doi:10.5061/dryad.mb0sd
# These functions are for use with the "variable rates" trees.
#
# Each simulated tree in MEA is associated with the following:
# 	An R data object, e.g., sim_4.Rda
#	The actual tree, e.g., sim_4.tre
# When you load the .Rda data file, you will notice an object "pruned_tree"
# 	appear in your workspace. This object is the phylogeny with an auxiliary
#   component "full_process", which contains the rate parameters for each 
#   branch. The functions below process this matrix to generate the mean branch 
#   rates, or the tip rates. 
#
#   Argument x to each function is the pruned_tree object
#
#	Return values: 
# 		branchRatesMEA: matrix of mean rates for each branch, with row 
#       names equal to the ape format node numbers for the tree stored 
#       in the pruned_tree object
#
#		tipRatesMEA: matrix of tip rates, with row 
#       names equal to the ape format node numbers for the tree stored 
#       in the pruned_tree object. Also includes the "regime_index", the 
#       index value of the rate regime to which each tip belongs


branchRatesMEA <- function(x){

	fx <- function(z){
		zz <- unlist(z)
		return(as.numeric(z[length(z)]))
	}
	
	zz <- x$full_process
	
	regime <- unlist(lapply(zz$states, fx))
	
	lambda <- numeric(nrow(zz))
	mu <- numeric(nrow(zz))
	
	for (i in 1:nrow(zz)){
 
	
		tmpT <- unlist(zz$transition_times[[i]])
		tmpL <- unlist(zz$speciation_rates[[i]])
		tmpM <- unlist(zz$extinction_rates[[i]])
		
		if (length(tmpT) != length(tmpL) | length(tmpT) != length(tmpM)){
			stop("full_process matrix decomp error")
		}
		
		if (length(tmpT) == 1){
			lambda[i] <- tmpL
			mu[i] <- tmpM
		}else{
			bl <- zz$end_time[i] - zz$start_time[i]
 			splits <- c(zz$start_time[i], tmpT[2:length(tmpT)], zz$end_time[i])
			wts <- diff(splits) / bl
	 
			for (k in 1:length(wts)){
				lambda[i] <- lambda[i] + tmpL[k]*wts[k]
				mu[i] <- mu[i] + tmpM[k]*wts[k]
			}
			
		}
		
	}
	
	mm <- cbind(lambda, mu)
	colnames(mm) <- c("lambda_branchmean", "mu_branchmean")

	rownames(mm) <- x$edge[,2]
	return(mm)
}


tipRatesMEA <- function(x){
	
	fx <- function(z){
		zz <- unlist(z)
		return(as.numeric(z[length(z)]))
	}
	
	regime <- unlist(lapply(x$full_process$states, fx))
	names(regime) <- x$edge[,2]
	regime <- regime[as.character(1:length(x$tip.label))]
	
	lrates <- x$full_process$current_speciation_rate
	names(lrates) <- x$edge[,2]
	lrates <- lrates[as.character(1:length(x$tip.label))]
	
	mrates <- x$full_process$current_extinction_rate
	names(mrates) <- x$edge[,2]
	mrates <- mrates[as.character(1:length(x$tip.label))]
 
	mat <- cbind(regime, lrates, mrates)
	rownames(mat) <- names(lrates)
	colnames(mat) <- c("regime_index", "lambda_tip", "mu_tip")
	
 	return(mat)
}
