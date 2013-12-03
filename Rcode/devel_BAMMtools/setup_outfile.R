setup_outfile <- function(fname, parnames, parvalues, comments){
	cat('#\n#\tFile generated ', as.character(Sys.time()), '\n', file = fname, append=F);
	cat('#\t\t\tby BAMMtools::autotune_BAMM_MCMC\n', file = fname, append=T);	
	cat('#\n#\tThis file contains autotuned parameters\n', file = fname, append=T);
	cat('#\tfor a BAMM analysis\n#\n', file = fname, append=T);
	cat('#\tCopy the block of parameter below\n', file = fname, append=T);	
	cat('#\tinto the MCMC_OPERATORS section of your controlfile\n', file = fname, append=T);		
	cat('#\tEnsure that the default MCMC_SCALING_OPERATORS section is deleted\n', file = fname, append=T);	
	cat('#\tstart MCMC_SCALING_OPERATORS block <<<<\n\n', file = fname, append=T);	
	for (i in 1:length(parnames)){
		cat(parnames[i], ' = ', parvalues[i], sep='', file = fname, append=T);
		if (comments[i] != ''){
			cat('   #', comments[i], sep='', file = fname, append=T);				
		}
		cat('\n', file = fname, append=T);		
	}
	cat('\n#\t>>>>>end MCMC_SCALING_OPERATORS block', file = fname, append=T);	
	cat('\n#\tNote: if you received a warning message after a parameter,', file = fname, append=T);	
	cat('\n#\tthis just means that the auto-tuning could not identify a best target value', file = fname, append=T);	
	cat('\n#\tfor this scaling operator. You do not need to worry about this.', file = fname, append=T);	
	cat('\n#\tHowever, if MCMC is performing poorly, you may consider manually', file = fname, append=T);	
	cat('\n#\tIncreasing or decreasing the value of the parameter', file = fname, append=T);	
}
