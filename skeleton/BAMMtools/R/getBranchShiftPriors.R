
getBranchShiftPriors <- function(phy, priordata){
	
	if (class(priordata) == 'character'){
		prior <- read.csv(priordata, header=T);
	}else if (class(priordata) == 'data.frame'){
		prior <- priordata;
	}else{
		stop("invalid priordata argument (wrong class) in getBranchShiftPriors\n");
	}
	
	tx <- cumsum(table(priordata$N_shifts) / nrow(priordata));
	tx <- as.numeric(names(tx[tx >=  0.95][1]));
	
	wts <- phy$edge.length / sum(phy$edge.length);

	obj <- phy;	
	obj$edge.length <- wts * tx;
	obj$criterion <- 0.95;
	class(obj) <- 'branchprior';
	return(obj);
}

