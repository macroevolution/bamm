#############################################################
#
#	addBAMMlegend(x, corners = c(0,1,0,10), side = 3, nTicks = 2, axisOffset = 0.4, ...)
#
#		x = saved plot.bammdata object
#		corners = coordinates for legend c(xmin, xmax, ymin, ymax)
#		side = side for tick marks, see axis() documentation
#		nTicks = number of ticks, outside of min and max
#		axisOffset = distance from color bar for labels
#		... additional parameters to be passed to axis()
#



addBAMMlegend <- function(x, corners = c(0,1,0,10), side = 3, nTicks = 2, direction = 'auto', axisOffset = 0.4, ...) {
	#corners xmin,xmax,ymin,ymax
	
	if (!identical(names(x), c('coords', 'colorbreaks', 'palette', 'colordens'))) {
		stop("x must be a saved plot.bammdata object.")
	}
	
	if (!direction %in% c('auto', 'vertical', 'horizontal')) {
		stop("direction must be auto, vertical or horizontal.");
	}
	
	colorbreaks <- x$colorbreaks
	pal <- x$palette

	n <- length(colorbreaks);
	
	if (direction == 'auto') {
		if ((corners[2] - corners[1]) >= (corners[4] - corners[3])) {
			direction <- 'horizontal';
		} else {
			direction <- 'vertical';
		}
	}

	if (direction == 'horizontal') {
		x <- seq(from = corners[1], to = corners[2], length.out = n);
		width <- corners[3:4];
	} else {
		x <- seq(from = corners[3], to = corners[4], length.out = n);
		width <- corners[1:2];
	}
	
	#get bin coordinates
	x <- rep(x,each = 2);
	x <- x[-c(1,length(x))];
	x <- matrix(x, ncol = 2, byrow = TRUE);
	
	#find tick locations
	#get equivalent rate bins
	z <- rep(colorbreaks,each = 2);
	z <- z[-c(1,length(z))];
	z <- matrix(z, ncol = 2, byrow = TRUE);

	tx <- trunc(seq(from = 1, to = nrow(x), length.out = nTicks + 2));
	tickLocs <- x[tx,1]
	tx <- z[tx,1]
	tickLocs[length(tickLocs)] <- max(x[,2])
	tx[length(tx)] <- max(z[,2])	
	
	#plot bar
	if (direction == 'horizontal') {
		rect(xleft = x[,1], ybottom = width[1], xright = x[,2], ytop = width[2], border = pal, col = pal);
	} else {
		rect(xleft = width[1], ybottom = x[,1], xright = width[2], ytop = x[,2], border = pal, col = pal);
	}
	
	#add tickmarks
	if (side == 1) { #bottom
		axis(side, at = tickLocs, pos = corners[3] - axisOffset, labels = signif(tx, 2), las = 1, ...);
	} 
	if (side == 3) { #top
		axis(side, at = tickLocs, pos = corners[4] + axisOffset, labels = signif(tx, 2), ...);
	}
	if (side == 2) { #left
		axis(side, at = tickLocs, pos = corners[1] - axisOffset, labels = signif(tx, 2), ...);
	}
	if (side == 4) { #right
		axis(side, at = tickLocs, pos = corners[2] + axisOffset, labels = signif(tx, 2), ...);
	}
}
	
