############################################
#	Internal function called by plot.bammdata(...)
#
#
colorMap = function(x, pal, NCOLORS, show=FALSE, ...) {
	dpal = c('BrBG','PRGn','PiYG','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral');
	if (length(pal) == 3) {
		colpalette = colorRampPalette(pal,space='Lab')(NCOLORS);	
	}
	else if (pal %in% dpal) {
		colpalette = colorRampPalette(rev(brewer.pal(3,pal)),space='Lab')(NCOLORS);
	}
	else if (tolower(pal) == "set1") {
		colpalette = colorRampPalette(c("darkgreen","yellow2","red"),space='Lab')(NCOLORS);
	}
	else if (tolower(pal) == "set2") {
		colpalette = colorRampPalette(c("darkgreen","pink","magenta4"),space='Lab')(NCOLORS);
	}
	else if (pal == 'temperature') {
		colpalette = richColors(NCOLORS);	
	}
	kde = density(x);
	bks = quantile(x, seq(0,1,length.out=(NCOLORS+1)));
	colset = numeric(length(x));
	coldens = numeric(length(kde$x));
	for (i in 2:length(bks)) {
        if (i == 2) {
            colset[x < bks[2]] = colpalette[1];
            coldens[kde$x < bks[2]] = colpalette[1];
        }
        else if (i == length(bks)) {
            colset[x >= bks[length(bks)-1]] = colpalette[length(bks)-1];
            coldens[kde$x >= bks[length(bks)-1]] = colpalette[length(bks)-1];
        }
        else {
            colset[x >= bks[i-1] & x < bks[i]] = colpalette[i-1];
        	coldens[kde$x >= bks[i-1] & kde$x < bks[i]] = colpalette[i-1];
        }
    }
	
	bks = rep(bks, each=2);
	bks = bks[-c(1,length(bks))];
	bks = matrix(bks, ncol=2, byrow=TRUE);
	bks = data.frame(colpalette,bks,stringsAsFactors=FALSE);
	colnames(bks) = c("Color","Lb","Hb");
		
	coldens = data.frame(kde$x,kde$y,coldens,stringsAsFactors=FALSE);
	if (show) {
		plot(coldens[,1],coldens[,2],type="n")
		for(i in 1:nrow(coldens)) {
			lines(rep(coldens[i,1],2),c(0,coldens[i,2]), col=coldens[i,3],...);	
		}
	}
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