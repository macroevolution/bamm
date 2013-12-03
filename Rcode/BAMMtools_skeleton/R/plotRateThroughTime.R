plotRateThroughTime <-
function(ephy, useMedian = F, intervals=seq(from = 0,to = 1,by = 0.01), ratetype = 'speciation', nBins = 100, smooth = F, smoothParam = 0.20, opacity = 0.01, intervalCol='blue', avgCol='red',start.time = NULL, end.time = NULL, node = NULL, nodetype='include', plot = T, cex.axis=1, cex=1.3, xline=3.5, yline=3.5, mar=c(6,6,1,1), xticks=5, yticks=5, xlim='auto', ylim='auto',add=F){
	
	if (!'bamm-data' %in% class(ephy) & !'bamm-ratematrix' %in% class(ephy)){
		stop("ERROR: Object ephy must be of class bamm-data\n or bamm-ratematrix.");
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

	if ('bamm-data' %in% class(ephy)){
		#get rates through binned time
		rmat <- getRateThroughTimeMatrix(ephy, start.time = start.time, end.time = end.time,node = node, nslices = nBins, nodetype=nodetype);
	}
	if ('bamm-ratematrix' %in% class(ephy)){
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
	if (useMedian == F){
		avg <- colMeans(rate);
	}
	if (useMedian == T){
		avg <- unlist(apply(rate,2,median));
	}
	
	#apply loess smoothing to intervals
	if (smooth == T){
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
	if (plot == T){
		if (add == F){
			plot.new();
			par(mar=mar);
			if (xlim == 'auto' & ylim == 'auto'){
				plot.window(xlim=c(maxTime, 0), ylim=c(0 , max(poly[[1]][,2])));
				xMin <- 0;
				xMax <- maxTime;
				yMin <- 0;
				yMax <- max(poly[[1]][,2]);
			}
			if (xlim != 'auto'){
				plot.window(xlim = xlim, ylim=c(0 , max(poly[[1]][,2])));
				xMin <- xlim[2];
				xMax <- xlim[1];
				yMin <- 0;
				yMax <- max(poly[[1]][,2]);
			}
			if (ylim != 'auto'){
				plot.window(xlim=c(maxTime, 0), ylim=ylim);
				xMin <- 0;
				xMax <- maxTime;
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

		axis(at=c(-1,round(seq(xMin, 1.3*xMax, length.out=xticks+1))), labels = c(-1,round(seq(xMin, 1.3*xMax, length.out=xticks+1))), cex.axis = cex.axis, side = 1);
		axis(at=c(-0.2,seq(yMin, 1.2*yMax, length.out=yticks+1)), labels = c(-0.2,round(seq(yMin, 1.2*yMax, length.out=yticks+1),digits=1)), las=1, cex.axis = cex.axis, side = 2);

		mtext(side = 1, text = 'Time since present', line = xline, cex = cex);
		mtext(side = 2, text = ratelabel, line = yline, cex = cex);
	}
	if (plot == F){
		return(list(poly = poly,avg = avg,times = rmat$time));
	}
}
