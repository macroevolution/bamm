############################################
#	Internal function called by plot.bammdata(...)
#
#
colorMap <- function(x, pal, breaks) {
    dpal <- list(
BrBG=rev(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30")),
PiYG=rev(c("#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef","#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221","#276419")),
PRGn=rev(c("#40004b","#762a83","#9970ab","#c2a5cf","#e7d4e8","#f7f7f7","#d9f0d3","#a6dba0","#5aae61","#1b7837","#00441b")),
PuOr=rev(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6","#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788","#2d004b")),
RdBu=rev(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061")),
RdYlBu=rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")),
BuOr=c("#002bff","#1a66ff","#3399ff","#66CCff","#99eeff","#ccffff","#ffffcc","#ffee99","#ffee66","#ff9933","#ff661a","#ff2b00"),
BuOrRd=c("#085aff","#3377ff","#5991ff","#8cb2ff","#bfd4FF","#e6eeff","#f7faff","#ffffcc","#ffff99","#ffff00","#ffcc00","#ff9900","#ff6600","#ff0000"),
DkRdBu=c("#2a0bd9","#264eff","#40a1ff","#73daff","#abf8ff","#e0ffff","#ffffbf","#ffe099","#ffad73","#f76e5e","#d92632","#a60021"),
BuDkOr=c("#1f8f99","#52c4cc","#99faff","#b2fcff","#ccfeff","#e6ffff","#ffe6cc","#ffca99","#ffad66","#ff8f33","#cc5800","#994000"),
GnPu=c("#005100","#008600","#00bc00","#00f100","#51ff51","#86ff86","#bcffbc","#ffffff","#fff1ff","#ffbcff","#ff86ff","#ff51ff","#f100f1","#bc00bc","#860086","#510051"),
RdYlGn=rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850","#006837")),
Spectral=c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2")
);
    	
	NCOLORS <- length(breaks)-1;
	if (length(pal) >= 3) {
		colpalette <- colorRampPalette(pal,space='Lab')(NCOLORS);	
	}
	else if (pal %in% names(dpal)) {
		colpalette <- colorRampPalette(dpal[[pal]],space='Lab')(NCOLORS);
	}
	else if (tolower(pal) == "temperature") {
		colpalette <- richColors(NCOLORS);	
	}
	else if (tolower(pal) == "terrain") {
		colpalette <- terrain.colors(NCOLORS);
	}
	else {
		stop("Unrecognized color palette specification");
	}
	kde <- density(x, from=min(x), to=max(x));
	colset <- numeric(length(x));
	coldens <- numeric(length(kde$x));
	for (i in 2:length(breaks)) {
        if (i == 2) {
            colset[x < breaks[2]] <- colpalette[1];
            coldens[kde$x < breaks[2]] <- colpalette[1];
        }
        else if (i == length(breaks)) {
            colset[x >= breaks[length(breaks)-1]] <- colpalette[length(breaks)-1];
            coldens[kde$x >= breaks[length(breaks)-1]] <- colpalette[length(breaks)-1];
        }
        else {
            colset[x >= breaks[i-1] & x < breaks[i]] <- colpalette[i-1];
        	coldens[kde$x >= breaks[i-1] & kde$x < breaks[i]] <- colpalette[i-1];
        }
    }
	coldens <- data.frame(kde$x,kde$y,coldens,stringsAsFactors=FALSE);
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