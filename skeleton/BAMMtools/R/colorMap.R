############################################
#	Internal function called by plot.bammdata(...)
#
#
colorMap = function(x, pal, breaks) {
	dpal = c('BrBG','PRGn','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral');
	NCOLORS = length(breaks)-1;
	if (length(pal) == 3) {
		colpalette = colorRampPalette(pal,space='Lab')(NCOLORS);	
	}
	else if (pal %in% dpal) {
		if (require(RColorBrewer) == FALSE) {
			stop("You specified a palette the requires the package RColorBrewer, which cannot be loaded");
		}
		colpalette = colorRampPalette(rev(brewer.pal(3,pal)),space='Lab')(NCOLORS);
	}
	else if (tolower(pal) == "set1") {
		colpalette = colorRampPalette(c("darkgreen","yellow2","red"),space='Lab')(NCOLORS);
	}
	else if (tolower(pal) == "set2") {
		colpalette = colorRampPalette(c("darkgreen","pink","magenta4"),space='Lab')(NCOLORS);
	}
	else if (tolower(pal) == "temperature") {
		colpalette = richColors(NCOLORS);	
	}
	else if (tolower(pal) == "terrain") {
		colpalette = terrain.colors(NCOLORS);
	}
	else {
		stop("Unrecognized color palette specification");
	}
	kde = density(x, from=min(x), to=max(x));
	colset = numeric(length(x));
	coldens = numeric(length(kde$x));
	for (i in 2:length(breaks)) {
        if (i == 2) {
            colset[x < breaks[2]] = colpalette[1];
            coldens[kde$x < breaks[2]] = colpalette[1];
        }
        else if (i == length(breaks)) {
            colset[x >= breaks[length(breaks)-1]] = colpalette[length(breaks)-1];
            coldens[kde$x >= breaks[length(breaks)-1]] = colpalette[length(breaks)-1];
        }
        else {
            colset[x >= breaks[i-1] & x < breaks[i]] = colpalette[i-1];
        	coldens[kde$x >= breaks[i-1] & kde$x < breaks[i]] = colpalette[i-1];
        }
    }
	coldens = data.frame(kde$x,kde$y,coldens,stringsAsFactors=FALSE);
	return(list(cols = colset, colsdensity = coldens));
}

# colorMap = function(x, pal, NCOLORS)
# {
	# dpal = c('BrBG','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral');
	# colset = numeric(length(x));
	# if(length(pal) == 3)
	# {
		# colpalette = colorRampPalette(pal,space='Lab')(NCOLORS);	
	# }
	# else if(pal %in% dpal)
	# {
		# colpalette = colorRampPalette(rev(brewer.pal(3,pal)),space='Lab')(NCOLORS);
	# }
	# else if(pal == 'temperature')
	# {
		# colpalette = richColors(NCOLORS);	
	# }
	
	# bks = quantile(x, seq(0,1,length.out=(NCOLORS+1)));
	# for(i in 2:length(bks))
	# {
		# if(i == 2)
		# {
			# colset[x < bks[2]] = colpalette[1];	
		# }
		# else if(i == length(bks))
		# {
			# colset[x >= bks[length(bks)-1]] = colpalette[length(bks)-1];
		# }
		# else
		# {
			# colset[x >= bks[i-1] & x < bks[i]] = colpalette[i-1];
		# }
	# }
	# return(colset);
# }