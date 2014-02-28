
summary.credibleshiftset <- function(x){

	conf <- round(100*x$set.limit);
	cat('\n', conf, '% credible set of rate shift configurations sampled with BAMM');
	cat('\n\nDistinct shift configurations in credible set: ', x$number.distinct);
	
	mm <- min(c(9, x$number.distinct));
	if (x$number.distinct > 9){
		omitted <- x$number.distinct - 9;
	}
	
	xvec <- numeric(x$number.distinct);
	dd <- data.frame(rank=1:x$number.distinct, probability = xvec, cumulative=xvec, N_shifts=xvec);
	for (i in 1:nrow(dd)){
		dd$probability[i] <- x$frequency[i];
		dd$cumulative[i] <- x$cumulative[i];
		dd$N_shifts[i] <- length(x$shiftnodes[[i]]);
	}
	
	
	cat('\n\nFrequency of', mm, 'shift configurations with highest posterior probability:\n');
	cat('\n\tRank\t\tProb\tCumulative\t\tCoreShifts\n')
	for (i in 1:mm){
		cat('\t', format(i, justify='left'), '\t\t', format(round(x$frequency[i], 3), nsmall=3), '\t');
		cat(round(x$cumulative[i], 3), '\t\t\t', length(x$shiftnodes[[i]]), '\n');
	}	
	
	if (x$number.distinct > 9){
		cat('\n...omitted', omitted, 'additional distinct shift configurations\n');
		cat('from the credible set. You can access the full set from your \n');
		cat('credibleshiftset object\n');
	}
	invisible(dd);
}


