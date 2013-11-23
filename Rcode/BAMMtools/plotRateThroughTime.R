
#############################################################
#
# plotRateThroughTime <- function(...)
#
# useMedian = boolean, will plot median if TRUE, mean if FALSE.
# intervals if NULL, no intervals will be plotted, otherwise a vector of quantiles must be supplied (these will define shaded polygons)
# ratetype = 'speciation' or 'extinction' or 'netdiv' or 'trait'
# nBins = number of time slices used to generate rates through time
# smooth = boolean whether or not to apply loess smoothing
# smoothParam = loess smoothing parameter, ignored if smooth = F
# opacity = opacity of color for interval polygons
# intervalCol = transparent color for interval polygons
# avgCol = color for mean/median line
# start.time = start time to be fed to getRateThroughTimeMatrix
# end.time = end time to be fed to getRateThroughTimeMatrix
# node = if supplied, the clade descended from this node will be used.
# nodetype = supplied to getRateThroughTimeMatrix
# plot = boolean: if TRUE, a plot will be returned, if FALSE, the data for the plot will be returned. 
#
plotRateThroughTime <- function(ephy, useMedian = F, intervals=seq(from = 0,to = 1,by = 0.01), ratetype = 'speciation', nBins = 100, smooth = F, smoothParam = 0.20, opacity = 0.01, intervalCol='blue', avgCol='red',start.time = NULL, end.time = NULL, node = NULL, nodetype='include', plot = T){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	if (!is.logical(useMedian)){
		stop('ERROR: useMedian must be either TRUE or FALSE.');
	}
	if (class(intervals)!='numeric' & class(intervals)!='NULL'){
		stop("ERROR: intervals must be either 'NULL' or a vector of quantiles.");
	}

	#get rates through binned time
	rmat <- getRateThroughTimeMatrix(ephy, start.time = start.time, end.time = end.time,node = node, nslices = nBins, nodetype=nodetype);

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

		poly<-list();
		q1<-1;
		q2<-nrow(mm);
		repeat{
			if (q1 >= q2) {break}
			a<-as.data.frame(cbind(rmat$times,mm[q1,]));
			b<-as.data.frame(cbind(rmat$times,mm[q2,]));
			b<-b[rev(rownames(b)),];
			colnames(a)<-colnames(b)<-c('x','y');
			poly[[q1]]<-rbind(a,b);
			q1<-q1+1;
			q2<-q2-1;
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
			poly[[i]][1:nrow(poly[[i]])/2,2] <- loess(poly[[i]][1:nrow(poly[[i]])/2,2] ~ poly[[i]][1:nrow(poly[[i]])/2,1],span = smoothParam)$fitted;
			poly[[i]][(nrow(poly[[i]])/2):nrow(poly[[i]]),2] <- loess(poly[[i]][(nrow(poly[[i]])/2):nrow(poly[[i]]),2] ~ poly[[i]][(nrow(poly[[i]])/2):nrow(poly[[i]]),1],span = smoothParam)$fitted;
		}
		avg <- loess(avg ~ rmat$time,span = smoothParam)$fitted;
	}

	#begin plotting
	if (plot == T){
		plot.new();
		plot.window(xlim=c(maxTime, 0), ylim=c(0 , max(poly[[1]][,2])));
	
		#plot intervals
		if (!is.null(intervals)){
			for (i in 1:length(poly)){
				polygon(x=maxTime - poly[[i]][,1],y=poly[[i]][,2],col=transparentColor(intervalCol,opacity),border=NA);
			}
		}
		lines(x = maxTime - rmat$time, y = avg, lwd = 3, col = avgCol);

		axis(at=seq(0, maxTime + 0.3*maxTime, by = 5), cex.axis = 1, side = 1);
		axis(at=seq(-0.2, max(rate) + 0.2*max(rate), by=0.1), las=1, cex.axis = 1, side = 2);
		mtext(side = 1, text = 'Time since present', line = 3, cex = 1.1);
		mtext(side = 2, text = ratelabel, line = 3, cex = 1.1);
	}
	if (plot == F){
		return(list(poly = poly,avg = avg,times = rmat$time))
	}
}
