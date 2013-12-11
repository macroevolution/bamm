
#############################################################
#
#	plotRateThroughTime <- function(...)
#
#	ephy = object of class bamm-data or bamm-ratematrix
#		if bamm-ratematrix, start.time, end.time, node, nslices, nodetype are not used.
#	useMedian = boolean, will plot median if TRUE, mean if FALSE.
#	intervals if NULL, no intervals will be plotted, otherwise a vector of quantiles must be supplied (these will define shaded polygons)
#	ratetype = 'speciation' or 'extinction' or 'netdiv' or 'trait'
#	nBins = number of time slices used to generate rates through time
#	smooth = boolean whether or not to apply loess smoothing
#	smoothParam = loess smoothing parameter, ignored if smooth = F
#	opacity = opacity of color for interval polygons
#	intervalCol = transparent color for interval polygons
#	avgCol = color for mean/median line
#	start.time = start time to be fed to getRateThroughTimeMatrix
#	end.time = end time to be fed to getRateThroughTimeMatrix
#	node = if supplied, the clade descended from this node will be used.
#	nodetype = supplied to getRateThroughTimeMatrix
#	plot = boolean: if TRUE, a plot will be returned, if FALSE, the data for the plot will be returned. 
#	xticks = number of ticks on the x-axis.
#	yticks = number of ticks on the y-axis.
#	xlim = vector of length 2 with min and max times for x axis. X axis is time since present, so if plotting till the present, xlim[2]==0. Can also be 'auto'.
#	ylim = vector of length 2 with min and max rates for y axis. Can also be 'auto'. 
#	add = boolean: should rates be added to an existing plot
#
#	+ several undocumented args to set plot parameters: mar, cex, xline, yline, etc.
#	

plotRateThroughTime <- function(ephy, useMedian = FALSE, intervals=seq(from = 0,to = 1,by = 0.01), ratetype = 'speciation', nBins = 100, smooth = FALSE, smoothParam = 0.20, opacity = 0.01, intervalCol='blue', avgCol='red',start.time = NULL, end.time = NULL, node = NULL, nodetype='include', plot = TRUE, cex.axis=1, cex=1.3, xline=3.5, yline=3.5, mar=c(6,6,1,1), xticks=5, yticks=5, xlim='auto', ylim='auto',add=FALSE){
	
	if ('bammdata' != class(ephy) & 'bamm-ratematrix' != class(ephy)){
		stop("ERROR: Object ephy must be of class bammdata or bamm-ratematrix.\n");
	}
	if (!is.logical(useMedian)){
		stop('ERROR: useMedian must be either TRUE or FALSE.');
	}
	if (class(intervals)!='numeric' & class(intervals)!='NULL'){
		stop("ERROR: intervals must be either 'NULL' or a vector of quantiles.");
	}
	if (!ratetype %in% c('speciation','trait','extinction','netdiv')){
		stop("ERROR: ratetype must be either 'speciation', 'extinction','netdiv' or 'trait'.");
	}
	if (!is.logical(smooth)){
		stop('ERROR: smooth must be either TRUE or FALSE.');
	}

	if ('bammdata' == class(ephy)){
		#get rates through binned time
		rmat <- getRateThroughTimeMatrix(ephy, start.time = start.time, end.time = end.time,node = node, nslices = nBins, nodetype=nodetype);
	}
	if ('bamm-ratematrix' == class(ephy)){
		#use existing rate matrix
		rmat <- ephy;
	}

	#set appropriate rates
	if (ratetype != 'speciation' & ratetype != 'extinction' & ratetype != 'netdiv' & ratetype != 'trait'){
		stop("ERROR: ratetype must be 'speciation', 'extinction', 'netdiv' or 'trait'.\n");
	}
	if (ratetype == 'speciation'){
		rate <- rmat$lambda;
		ratelabel <- 'Speciation';
	}
	if (ratetype == 'trait'){
		rate <- rmat$lambda;
		ratelabel <- 'BM rate';
	}
	if (ratetype == 'extinction'){
		rate <- rmat$mu;
		ratelabel <- 'Extinction';
	}
	if (ratetype == 'netdiv'){
		rate <- rmat$lambda - rmat$mu;
		ratelabel <- 'Net diversification';
	}

	#generate coordinates for polygons
	maxTime <- max(rmat$times);
	if (!is.null(intervals)){
		mm <- apply(rate, MARGIN = 2, quantile, intervals);

		poly <- list();
		q1 <- 1;
		q2 <- nrow(mm);
		repeat{
			if (q1 >= q2) {break}
			a <- as.data.frame(cbind(rmat$times,mm[q1,]));
			b <- as.data.frame(cbind(rmat$times,mm[q2,]));
			b <- b[rev(rownames(b)),];
			colnames(a) <- colnames(b) <- c('x','y');
			poly[[q1]] <- rbind(a,b);
			q1 <- q1 + 1;
			q2 <- q2 - 1;
		}
	}

	#Calculate averaged data line
	if (useMedian == FALSE){
		avg <- colMeans(rate);
	}
	if (useMedian == TRUE){
		avg <- unlist(apply(rate,2,median));
	}
	
	#apply loess smoothing to intervals
	if (smooth == TRUE){
		for (i in 1:length(poly)){
			p <- poly[[i]];
			rows <- nrow(p);
			p[1:rows/2,2] <- loess(p[1:rows/2,2] ~ p[1:rows/2,1],span = smoothParam)$fitted;
			p[(rows/2):rows,2] <- loess(p[(rows/2):rows,2] ~ p[(rows/2):rows,1],span = smoothParam)$fitted;
			poly[[i]] <- p;
		}
		avg <- loess(avg ~ rmat$time,span = smoothParam)$fitted;
	}

	#begin plotting
	if (plot == TRUE){
		if (add == FALSE){
			plot.new();
			par(mar=mar);
			if (unique(xlim == 'auto') & unique(ylim == 'auto')){
				plot.window(xlim=c(maxTime, 0), ylim=c(0 , max(poly[[1]][,2])));
				xMin <- maxTime;
				xMax <- 0;
				yMin <- 0;
				yMax <- max(poly[[1]][,2]);
			}
			if (unique(xlim != 'auto') & unique(ylim == 'auto')){
				plot.window(xlim = xlim, ylim=c(0 , max(poly[[1]][,2])));
				xMin <- xlim[1];
				xMax <- xlim[2];
				yMin <- 0;
				yMax <- max(poly[[1]][,2]);
			}
			if (unique(xlim == 'auto') & unique(ylim != 'auto')){
				plot.window(xlim=c(maxTime, 0), ylim=ylim);
				xMin <- maxTime;
				xMax <- 0;
				yMin <- ylim[1];
				yMax <- ylim[2];
			}
		}
		#plot intervals
		if (!is.null(intervals)){
			for (i in 1:length(poly)){
				polygon(x=maxTime - poly[[i]][,1],y=poly[[i]][,2],col=transparentColor(intervalCol,opacity),border=NA);
			}
		}
		lines(x = maxTime - rmat$time, y = avg, lwd = 3, col = avgCol);

		axis(at=c(1.3*xMin,round(seq(xMin,xMax, length.out=xticks+1))), labels = c(1.3*xMin,round(seq(xMin, xMax, length.out=xticks+1))), cex.axis = cex.axis, side = 1);
		axis(at=c(-0.2,seq(yMin, 1.2*yMax, length.out=yticks+1)), labels = c(-0.2,round(seq(yMin, 1.2*yMax, length.out=yticks+1),digits=1)), las=1, cex.axis = cex.axis, side = 2);

		mtext(side = 1, text = 'Time since present', line = xline, cex = cex);
		mtext(side = 2, text = ratelabel, line = yline, cex = cex);
	}
	if (plot == FALSE){
		return(list(poly = poly,avg = avg,times = rmat$time));
	}
}
