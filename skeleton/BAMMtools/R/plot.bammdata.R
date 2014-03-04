plot.bammdata = function (x, method = "phylogram", vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, legend = FALSE, spex = "s", lwd = 1, cex = 1, pal = "RdYlBu", colorbreaks = NULL, par.reset = TRUE, ...) {
    if ("bammdata" %in% class(x)) {
    	if (attributes(x)$order != "cladewise") {
    		stop("Function requires tree in 'cladewise' order");
    	}
        phy <- as.phylo.bammdata(x);
    }
    else stop("Object ephy must be of class bammdata");
    if (!is.binary.tree(phy)) {
        stop("Function requires fully bifurcating tree");
    }
    if (any(phy$edge.length == 0)) {
        warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero");
    }
    if (!("dtrates" %in% names(x))) {
        x <- dtRates(x, 0.01);
    }
    if (is.null(colorbreaks)) {
    	colorbreaks <- assignColorBreaks(x$dtrates$rates, 64, spex);
    }
    if (x$type == "trait") {
        if (sum(is.na(x$dtrates$rates))) {
            warning(sprintf("Found %d NA phenotypic rates. Coercing to zero.", sum(is.na(x$dtrates$rates))));
            x$dtrates$rates[is.na(x$dtrates$rates)] <- 0;
        }
    	colorobj <- colorMap(x$dtrates$rates, pal, colorbreaks);
    }
    else if (x$type == "diversification") {
        if (sum(is.na(x$dtrates$rates[[1]]))) {
            warning(sprintf("Found %d NA speciation rates. Coercing to zero.", sum(is.na(x$dtrates$rates[[1]]))));
            x$dtrates$rates[[1]][is.na(x$dtrates$rates[[1]])] <- 0;
        }
        if (sum(is.na(x$dtrates$rates[[2]]))) {
            warning(sprintf("Found %d NA extinction rates. Coercing to zero.", sum(is.na(x$dtrates$rates[[2]]))));
            x$dtrates$rates[[2]][is.na(x$dtrates$rates[[2]])] <- 0;
        }
        if (tolower(spex) == "s") {
        	colorobj <- colorMap(x$dtrates$rates[[1]], pal, colorbreaks);
        }
        else if (tolower(spex) == "e") {
        	colorobj <- colorMap(x$dtrates$rates[[2]], pal, colorbreaks);
        }
        else {
        	colorobj <- colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]], pal, colorbreaks);
        }
    }
    else {
    	stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'");	
    }
    edge.color <- colorobj$cols;
    tH <- max(branching.times(phy));
    phy$begin <- x$begin;
    phy$end <- x$end;
    tau <- x$dtrates$tau;
    if (method == "polar") {
        ret <- setPolarTreeCoords(phy, vtheta, rbf);
        rb <- tH * rbf;
        p <- mkdtsegsPolar(ret$segs[-1,], tau, x$edge);
    }
    else if (method == "phylogram") {
        ret <- setPhyloTreeCoords(phy);
        p <- mkdtsegsPhylo(ret$segs[-1,], tau, x$edge);
    }
    else {
        stop("Unimplemented method");
    }
    x0 <- c(ret$segs[1,1], p[, 1]);
    x1 <- c(ret$segs[1,3], p[, 2]);
    y0 <- c(ret$segs[1,2], p[, 3]);
    y1 <- c(ret$segs[1,4], p[, 4]);
    offset <- table(p[, 5])[as.character(unique(p[, 5]))];
    arc.color <- c(edge.color[1], edge.color[match(unique(p[, 5]), p[, 5]) + offset]);
    edge.color <- c(edge.color[1], edge.color);
    if (show) {
        if (length(list(...))) {
            op <- par(no.readonly = TRUE);
            par(...);
        }
        plot.new();
        ofs <- 0;
        if (labels) {
            ofs <- max(nchar(phy$tip.label) * 0.03 * cex);
        }
        if (method == "polar") {
            plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs), ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, ofs), asp = 1);
            segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, lend = 2);
            arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + phy$end/tH), border = arc.color, lwd = lwd);
            if (labels) {
                for (k in 1:length(phy$tip.label)) {
                  text(ret$segs[-1, ][phy$edge[, 2] == k, 3],ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k],cex = cex, srt = (180/pi) * ret$arcs[-1,][phy$edge[, 2] == k, 1], adj = c(0, NA));
                }
            }
        }
        if (method == "phylogram") {
            plot.window(xlim = c(0, 1 + ofs), ylim = c(0, phy$Nnode * 1/(phy$Nnode + 1)));
            segments(x0[-1], y0[-1], x1[-1], y1[-1], col = edge.color[-1], lwd = lwd, lend = 2);
            isTip <- phy$edge[, 2] <= phy$Nnode + 1;
            isTip <- c(FALSE, isTip);
            segments(ret$arcs[!isTip, 1], ret$arcs[!isTip, 2], ret$arcs[!isTip, 3], ret$arcs[!isTip, 4], col = arc.color[!isTip], lwd = lwd, lend = 2);
            if (labels) {
                text(ret$segs[-1, ][phy$edge[, 2] <= phy$Nnode + 1, 3], ret$segs[-1, ][phy$edge[, 2] <= phy$Nnode + 1, 4], phy$tip.label, cex = cex, pos = 4, offset = 0.25);
            }
        }
        if (legend) {
            rateLegend(colorobj$colsdensity);
        }
    }
    index <- order(as.numeric(rownames(ret$segs)));
    if (method == "phylogram") {
        assign("last_plot.phylo", list(type = "phylogram", direction = "rightwards", Ntip = phy$Nnode + 1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 3], yy = ret$segs[index, 4], pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
    }
    else if (method == "polar") {
        assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode + 1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 3], yy = ret$segs[index, 4], theta = ret$segs[index, 5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
    }
    if (par.reset) {
        if (length(list(...))) {
            par(op);
        }
    }
    invisible(list(coords = ret$segs[-1, ], colorbreaks = colorbreaks, colordens = colorobj$colsdensity));
}

# plot.bammdata = function(ephy, method='phylogram', tau=0.01, index=NULL, vtheta=5, rbf=0.001, show=TRUE, labels=FALSE, multi=FALSE, hrates=FALSE, spex = "s", lwd=1, cex=1, ncolors=64, pal='temperature', ...) {
	# if ('bammdata' == class(ephy)) phy = as.phylo.bammdata(ephy) else stop('Object ephy must be of class bammdata\n');
	# if (!is.binary.tree(phy)) stop('Function requires fully bifurcating tree.');
	# if (any(phy$edge.length == 0)) warning('Tree contains zero length branches. Rates for these will be NA and coerced to zero');
	# if (!('dtrates' %in% names(ephy))) ephy = dtRates(ephy, tau, index);
	# if (ephy$type == "trait") {
		# if(sum(is.na(ephy$dtrates$rates))) {
			# warning(sprintf("Found %d NA phenotypic rates. Coercing to zero.",sum(is.na(ephy$dtrates$rates))));
			# ephy$dtrates$rates[is.na(ephy$dtrates$rates)] = 0;
		# }
		# edge.color = colorMap(ephy$dtrates$rates,pal,ncolors);
	# }
	# else if (ephy$type == "diversification") {
		# if(sum(is.na(ephy$dtrates$rates[[1]]))) {
			# warning(sprintf("Found %d NA speciation rates. Coercing to zero.",sum(is.na(ephy$dtrates$rates[[1]]))));
			# ephy$dtrates$rates[[1]][is.na(ephy$dtrates$rates[[1]])] = 0;
		# }
		# if(sum(is.na(ephy$dtrates$rates[[2]]))) {
			# warning(sprintf("Found %d NA extinction rates. Coercing to zero.",sum(is.na(ephy$dtrates$rates[[2]]))));
			# ephy$dtrates$rates[[2]][is.na(ephy$dtrates$rates[[2]])] = 0;
		# }
		# if (tolower(spex) == "s") {
			# edge.color = colorMap(ephy$dtrates$rates[[1]],pal,ncolors);
		# }
		# else if (tolower(spex) == "e") {
			# edge.color = colorMap(ephy$dtrates$rates[[2]],pal,ncolors);
		# }
		# else {
			# edge.color = colorMap(ephy$dtrates$rates[[1]] - ephy$dtrates$rates[[2]],pal,ncolors);
		# }
	# }
	# tH = max(branching.times(phy));
	# phy$begin = ephy$begin; phy$end = ephy$end;
	# if (method == 'polar') {
		# ret = setPolarTreeCoords(phy,vtheta,rbf);
		# rb = tH*rbf;
	# }	
	# else if (method == 'phylogram') {
		# ret = setPhyloTreeCoords(phy);
	# }
	# else {
		# stop('Unimplemented method');
	# }
	# x0 = ret$segs[,1];y0=ret$segs[,2];x1=ret$segs[,3];y1=ret$segs[,4];
	
	# tau = ephy$dtrates$tau;
	# p = cbind(x0[-1],y0[-1],x1[-1],y1[-1],phy$edge[,2]);
	# p = apply(p,1,mkdtsegs,tau,phy,tH);
	# p = do.call(rbind, p);
	# x0 = c(x0[1],p[,1]);x1=c(x1[1],p[,2]);y0=c(y0[1],p[,3]);y1=c(y1[1],p[,4]);
	# offset = table(p[,5])[as.character(unique(p[,5]))];
	# arc.color = c(edge.color[1],edge.color[match(unique(p[,5]),p[,5])+offset]);
	# edge.color = c(edge.color[1],edge.color);
	
	# if (show) {
		# if(length(list(...))) {
			# op = par(no.readonly=TRUE);
			# par(...);
		# }
		# plot.new(); ofs = 0;
		# if (labels) {
			# ofs = max(nchar(phy$tip.label) * 0.03 * cex);
		# }
		
		# if (method == 'polar')  {
			# plot.window(xlim=c(-1,1)+c(-rb,rb)+c(-ofs,ofs),ylim=c(-1,1)+c(-rb,rb)+c(-ofs,ofs),asp=1);
			# segments(x0,y0,x1,y1,col=edge.color,lwd=lwd,lend=2);	
			# arc(0,0,ret$arcs[,1],ret$arcs[,2],c(rb,rb+phy$end/tH),border=arc.color,lwd=lwd);
			# if(labels) {
				# for(k in 1:length(phy$tip.label)) {
					# text(ret$segs[-1,][phy$edge[,2]==k,3],ret$segs[-1,][phy$edge[,2]==k,4],phy$tip.label[k],cex=cex, srt = (180/pi)*ret$arcs[-1,][phy$edge[,2]==k,1],adj=c(0,NA));	
				# }
			# }
		# }
		# if (method == 'phylogram') {
			# plot.window(xlim=c(0,1+ofs),ylim=c(0,phy$Nnode*1/(phy$Nnode+1)));
			# segments(x0,y0,x1,y1,col=edge.color,lwd=lwd,lend=2);
			# isTip = phy$edge[,2] <= phy$Nnode+1; isTip = c(FALSE,isTip);
			# segments(ret$arcs[!isTip,1],ret$arcs[!isTip,2],ret$arcs[!isTip,3],ret$arcs[!isTip,4],col=arc.color[!isTip],lwd=lwd,lend=2);
			# if(labels) {
				# text(ret$segs[-1,][phy$edge[,2] <= phy$Nnode+1,3],ret$segs[-1,][phy$edge[,2] <= phy$Nnode+1,4], phy$tip.label, cex=cex, pos=4, offset = 0.25);
			# }
		# }
		# if(hrates) {
			# if (ephy$type == "trait") {
				# histRates(ephy$dtrates$rates,pal,ncolors);
			# }
			# else if (ephy$type == "diversification") {
				# if (tolower(spex) == "s") {
					# edge.color = histRates(ephy$dtrates$rates[[1]],pal,ncolors);
				# }
				# else if (tolower(spex) == "e") {
					# edge.color = histRates(ephy$dtrates$rates[[2]],pal,ncolors);
				# }
				# else {
					# edge.color = histRates(ephy$dtrates$rates[[1]] - ephy$dtrates$rates[[2]],pal,ncolors);
				# }	
			# }
		# }
	# }
	# index = order(as.numeric(rownames(ret$segs)));
	# if (method == 'phylogram') {
		# assign("last_plot.phylo", list(Ntip = phy$Nnode+1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,3], yy = ret$segs[index,4], pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
	# }
	# else if (method == 'polar') {
		# assign("last_plot.phylo", list(Ntip = phy$Nnode+1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index,3], yy = ret$segs[index,4], theta = ret$segs[index,5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv);
	# }
	# if(!multi) {
		# if(length(list(...))) par(op);
	# }
	# invisible(ret$segs[-1,]);
# }