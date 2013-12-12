summary.bammdata = function(ephy, opt="product", display=6)
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
	cat("Process number posterior distribution:\n\n");
	fev = sapply(ephy$eventData, nrow);
	disp = tabulate(fev); disp = disp/sum(disp);
	if (length(disp) < display)
	{	
		cat(signif(disp,2),"\n\n");
	}
	else
	{
		if(length(disp)-display == 1)
		{
			cat(disp[1:display],"... omitted 1 value\n\n");
		}
		else
		{
			cat(disp[1:display],"... omitted", length(disp)-display,"values\n\n");
		}
	}
	invisible(list(posterior = tabulate(fev), mscshiftnodes = shiftnodes));
}