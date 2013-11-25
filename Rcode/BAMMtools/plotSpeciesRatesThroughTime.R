#########################################################
#
#	plotSpeciesRatesThroughTime <- function(...)
#
#	Function to plot several species-specific rate trajectories together
#
#	start.time = time before present (defaults to full depth of tree)
#	useMedian = boolean, will plot median if TRUE, mean if FALSE
#	nbreaks = number of time slices to use
#	ratetype = 'speciation' or 'netdiv' or 'trait'
#	species = vector of species names to include
# 	intervals if NULL, no intervals will be plotted, otherwise a vector of quantiles must be supplied (these will define shaded polygons)
#	smooth = boolean, apply loess smoothing to curves
#	spCol = vector of colors that will be matched to the species vector (defaults to random color selection)
#	opacity = opacity level for plotting colored confidence intervals


plotSpeciesRatesThroughTime <- function(ephy, start.time=max(branching.times(ephy)), useMedian=F, nbreaks=10, ratetype='speciation', species, intervals=seq(from = 0,to = 1,by = 0.01), smooth=F, spCol=sample(colors(),2), opacity=0.01,smoothParam=0.20){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	if (!all(species %in% ephy$tip.label)){
		stop('Not all species names match those in ephy.');
	}
	
	if (ratetype=='speciation' | ratetype=='trait'){
		ndr <- F;
		if (ratetype == 'speciation'){ratelabel <- 'Speciation'}
		if (ratetype == 'trait'){ratelabel <- 'BM rate'}
	}
	if (ratetype == 'netdiv'){
		ndr <- T;
		ratelabel <- 'Net diversification';
	}
		
	bySp<-lapply(species,function(x) getSpeciesRateThroughTime(ephy, start.time = start.time, nbreaks = nbreaks, ndr = ndr, species = x,returnAll=T));
	names(bySp) <- species;
	
	#calculate average values
	if (useMedian == T){
		avgList<-lapply(bySp,function(x) apply(x,2,function(x) median(x)));
	}
	if (useMedian == F){
		avgList<-lapply(bySp,function(x) colMeans(x));
	}
	
	tend <- max(branching.times(ephy))*0.999;
	tstart <- max(branching.times(ephy)) - tend;
	tseq <- seq(tstart, tend, length.out = nbreaks);
	
	#calculate confidence intervals
	if (!is.null(intervals)){
		polySp<-list();
		for (i in 1:length(bySp)){
			mm <- apply(bySp[[i]], 2, quantile, intervals);
			poly <- list();
			q1 <- 1;
			q2 <- nrow(mm);
			repeat{
				if (q1 >= q2) {break}
				a <- as.data.frame(cbind(tseq,mm[q1,]));
				b <- as.data.frame(cbind(tseq,mm[q2,]));
				b <- b[rev(rownames(b)),];
				colnames(a) <- colnames(b) <- c('x','y');
				poly[[q1]] <- rbind(a,b);
				q1 <- q1 + 1;
				q2 <- q2 - 1;
			}
			polySp[[i]] <- poly;
		}
	}
	
	if (smooth == T){
		if (nbreaks < 30){
			cat('Insufficient breaks. Non-smoothed results returned.\n');
		}
		if (nbreaks >= 30){
			for (i in 1:length(polySp)){
				for (j in 1:length(polySp[[i]])){
					p <- polySp[[i]][[j]];
					rows <- nrow(p);
					p[1:rows/2,2] <- loess(p[1:rows/2,2] ~ p[1:rows/2,1],span = smoothParam)$fitted;
					p[(rows/2):rows,2] <- loess(p[(rows/2):rows,2] ~ p[(rows/2):rows,1],span = smoothParam)$fitted;
					polySp[[i]][[j]] <- p;
				}
			}
			for (i in 1:length(avgList)){
				avgList[[i]] <- loess(avgList[[i]] ~ tseq,span = smoothParam)$fitted;
			}
		}
	}
	
	maxRate <- max(unlist(lapply(polySp,function(x) lapply(x,function(x) x[,2]))));
	plot.new();
	plot.window(xlim=c(max(branching.times(ephy)), 0), ylim=c(0 , maxRate));

	if (!is.null(intervals)){
		for (i in 1:length(polySp)){
			for (j in 1:length(polySp[[i]])){
				polygon(x = tend - polySp[[i]][[j]][,1],y = polySp[[i]][[j]][,2], col = transparentColor(spCol[i],alpha = opacity), border = NA);
			}
			lines(x = tend - tseq,y = avgList[[i]], lwd=2, col = spCol[i]);
		}
	}
	if (is.null(intervals)){
		for (i in 1:length(avgList)){
			lines(x = tend - tseq,y = avgList[[i]],lwd = 2,col = spCol[i]);
		}
	}
	
	axis(at = seq(0, + 1.3*max(branching.times(ephy)), by = 5), cex.axis = 1, side = 1);
	axis(at = seq(-0.2, maxRate, by=0.2), las=1, cex.axis = 1, side = 2);
	mtext(side = 1, text = 'Time since present', line = 3, cex = 1.1);
	mtext(side = 2, text = ratelabel, line = 3, cex = 1.1);
}
