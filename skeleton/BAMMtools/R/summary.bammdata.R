summary.bammdata = function(object, display=10, ...) {
	cat("\nAnalyzed", length(object$eventData), "posterior samples\n");
#	shiftsindex <- maximumShiftCredibility(x,maximize=opt)$sampleindex;
#	shiftnodes <- getShiftNodesFromIndex(x, shiftsindex);
#	if (length(shiftnodes) > 1) {
#		cat("\nMaximum shift credibility tree has", length(shiftnodes), "shifts\n\n");
#		cat("Shifts occur leading to nodes:", shiftnodes,"\n\n");
#		cat("Optimality type:", opt,"\n\n");
#	}
#	else if (length(shiftnodes) == 1) {
#		cat("\nMaximum shift credibility tree has 1 shift\n\n");
#		cat("Shift occurs leading to node:", shiftnodes,"\n\n");
#		cat("Optimality type:", opt,"\n\n");
#	}
#	else {
#		cat("\nMaximum shift credibility tree has 0 shifts\n\n");
#		cat("Optimality type:", opt,"\n\n");
#	}
	cat("Shift posterior distribution:\n\n");
	fev <- sapply(object$eventData, nrow);
	disp <- tabulate(fev); disp = disp/sum(disp);
	disp <- data.frame(cbind(seq.int(0,length(disp)-1,1),signif(disp,2)));
	if (nrow(disp) <= display) {	
		write.table(format(disp,justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);	
	}
	else {
		if(nrow(disp)-display == 1) {
			write.table(format(disp[1:display,],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
			cat("... omitted 1 row\n\n");
		}
		else {
			wr <- which(disp[,2] == max(disp[,2]))[1];
			index <- c(max(1,floor(wr-display/2)), wr, min(ceiling(wr+display/2),nrow(disp)));
			if (index[1] > 1) {
				cat("... omitted",index[1]-1,"rows\n");
			}
			write.table(format(disp[index[1]:index[3],],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
			cat("... omitted", nrow(disp)-index[3]-1,"rows\n\n");
		}
	}
	cat("\nCompute credible set of shift configurations for more information:\n");
	cat("\tSee ?credibleShiftSet and ?getBestShiftConfiguration\n");
	
	xx <- table(fev);
	df <- data.frame(shifts = as.numeric(names(xx)), prob =  as.numeric(xx / (sum(xx))));
	
	invisible(df);
}
