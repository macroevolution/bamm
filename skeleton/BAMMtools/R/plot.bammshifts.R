## show.all.nodes: should plot core and "floater" nodes with different colors 
#	that can be specified by the user.
#    argument to specify layout , e.g., c(2,2) would specify a 2x2 panel?
#	 maybe print the number of shift configurations above the threshold that were 
#	 not plotted along with their frequency?
#	
#	add.freq.text = argument to add text to each plot indicating the frequency of the 
#		shift configuration

plot.bammshifts <- function(sc, ephy, plotmax=9, method='phylogram', pal = 'temperature', spex = "s", add.freq.text = TRUE, use.plot.bammdata = TRUE, send2pdf = FALSE){
	if (class(sc) != 'bammshifts') {
		stop('arg sc must be of class "bammshifts"');
	}
	if (class(ephy) != 'bammdata') {
		stop('arg ephy must be of class "bammdata"');
	}
	if (plotmax > 9 && send2pdf == FALSE) {
	    plotmax = 9;
	    cat("arg plotmax coerced to 9\n");
	}
	mm <- min(c(length(sc$frequency), plotmax));
	if (send2pdf) {
	    pdf("shiftconfig.pdf");
	}
	if (mm == 1) {
	    par(mfrow=c(1,1));
	} else if (mm <= 2) {
		par(mfrow=c(1,2));
	} else if (mm <= 4) {
		par(mfrow = c(2,2));
	} else if (mm <= 6) {
		par(mfrow=c(2,3));
	} else {
		par(mfrow=c(3,3));
	}
	cat("Omitted", max(length(sc$frequency),mm) - min(length(sc$frequency),mm), "plots\n");
	if (use.plot.bammdata) {
    	ephy = dtRates(ephy, 0.01);
	    colorbreaks = assignColorBreaks(ephy$dtrates$rates,spex=spex);
	}
	for (i in 1:mm) {
		tmp <- subsetEventData(ephy, index=sc$sampleset[[i]][1]);
		par(mar = c(2,2,2,2));
		if (use.plot.bammdata) {
    		plot.bammdata(tmp, method, pal=pal, colorbreaks=colorbreaks);
		}
		else {
		    if (method=="polar") method = "fan";
		    plot.phylo(as.phylo.bammdata(ephy),type=method,show.tip.label=FALSE);
		}
		if (add.freq.text) mtext(paste("f =",signif(sc$frequency[i],2)),3);
		box();
		cex = 2 + 8 * sc$marg.probs[as.character(getShiftNodesFromIndex(ephy,sc$sampleset[[i]][1]))];
		shiftnodes = getShiftNodesFromIndex(ephy,sc$sampleset[[i]][1]);
		shiftnode_parents = ephy$edge[which(ephy$edge[,2] %in% shiftnodes),1];
		
		isShiftNodeParent = integer(length(shiftnodes));		
		root = (shiftnode_parents == ephy$Nnode+2);
		isShiftNodeParent[root] = 1
		isShiftNodeParent[!root] = tmp$eventVectors[[1]][which(ephy$edge[,2] %in% shiftnode_parents[!root])];
		
		isShiftNode = which(tmp$eventData[[1]]$node %in% shiftnodes);
		
		AcDc = exponentialRate(tmp$eventData[[1]][isShiftNode, 2]-tmp$eventData[[1]][isShiftNodeParent, 2],tmp$eventData[[1]][isShiftNodeParent,3],tmp$eventData[[1]][isShiftNodeParent,4]) > tmp$eventData[[1]][isShiftNode, 3];
		
		bg = rep("blue", length(AcDc));
		bg[which(AcDc == FALSE)] = "red";
		
		#lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv);
        #XX = lastPP$xx[shiftnodes];
        #YY = lastPP$yy[shiftnodes];
        #points(XX, YY, pch = 21, cex = cex, col = 1, bg = transparentColor(bg,0.5)); 		
	    addBAMMshifts(tmp, method, 1, cex = cex, bg = transparentColor(bg,0.5), multi=TRUE);
	}
	if (send2pdf) dev.off();
}

