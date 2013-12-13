summary.bammdata = function(ephy, opt="product", display=10)
{
	cat("\nAnalyzed", length(ephy$eventData), "posterior samples\n");
	shiftsindex = maximumShiftCredibilityTree(ephy,maximize=opt)$sampleindex;
	shiftnodes = getShiftNodesFromIndex(ephy, shiftsindex);
	if (length(shiftnodes) > 1)
	{
		cat("\nMaximum shift credibility tree has", length(shiftnodes), "shifts\n\n");
		cat("Shifts occur leading to nodes:", shiftnodes,"\n\n");
		cat("Optimality type:", opt,"\n\n");
	}
	else if (length(shiftnodes) == 1)
	{
		cat("\nMaximum shift credibility tree has 1 shift\n\n");
		cat("Shift occurs leading to node:", shiftnodes,"\n\n");
		cat("Optimality type:", opt,"\n\n");
	}
	else
	{
		cat("\nMaximum shift credibility tree has 0 shifts\n\n");
		cat("Optimality type:", opt,"\n\n");
	}
	cat("Shift posterior distribution:\n\n");
	fev = sapply(ephy$eventData, nrow);
	disp = tabulate(fev); disp = disp/sum(disp);
	disp = data.frame(cbind(seq.int(0,length(disp)-1,1),signif(disp,2)));
	if (nrow(disp) <= display)
	{	
		write.table(format(disp,justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);	
	}
	else
	{
		if(nrow(disp)-display == 1)
		{
			write.table(format(disp[1:display,],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
			cat("... omitted 1 row\n\n");
		}
		else
		{
			write.table(format(disp[1:display,],justify="left",width=10), col.names=FALSE, row.names=FALSE, quote=FALSE);
			cat("... omitted", nrow(disp)-display,"rows\n\n");
		}
	}
	invisible(list(posterior = tabulate(fev), mscshiftnodes = shiftnodes));
}
