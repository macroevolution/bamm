library('colorspace');
library(gplots);

source('/Users/danrabosky/DanWork/bamm/devel/bamm/doc/sweave/BAMMtools.R');


getColorKey <- function(x, ncols = 16, units=2){
	
	if (ncols %% 2 != 0){
		stop("expecting even number of colors");
	}
	
	index <- (1:(ncols-1)) - (ncols/2);
	xv <- index * log(units);
	
	colset <- rich.colors(n=ncols);
	colmat <- matrix('', nrow=nrow(x), ncol=ncol(x));	
		
	for (i in 1:nrow(x)){
		for (j in 1:ncol(x)){
			if (!is.na(x[i,j])){
				if (x[i,j] <= xv[1]){
					colmat[i,j] <- colset[1];
				}else if (x[i,j] > xv[length(xv)]){
					colmat[i,j] <- colset[length(colset)];
				}else{
					for (k in 1:(length(xv)-1)){
						if (x[i, j] > xv[k] & x[i,j] <= xv[k+1]){
							colmat[i,j] <- colset[k];
						}
					}
				}

			}
		}
	}	
	
	return(colmat);
}

getColorBar <- function(ncols = 16, units=2){
	
	if (ncols %% 2 != 0){
		stop("expecting even number of colors");
	}
	
	index <- (1:(ncols-1)) - (ncols/2);
	xv <- index * log(units);
 
	colset <- rich.colors(n=ncols);
 	cx <- c(xv, max(xv)+log(units));
 
	return(data.frame(cuts=cx, cols=colset, stringsAsFactors=F));
}



cmat <- getColorKey(log(bfmat))

bfmat <- computeBayesFactors('post_mcmc_p50run1.txt', 'prior_mcmc_out.txt', modelset=0:100, constrain=F, burnin=0.25, threshold=0);





##########################
###### The figure!

quartz.options(height=7, width=10, dpi=72);
ll <- c(rep(1, 9), rep(2,3));
lmat <- matrix(ll, nrow=3, byrow=F);

plot.new();

layout(lmat);
plot.new();
par(mar=c(6,6,1,1));

plot.window(xlim=c(0,101), ylim=c(0,101), asp=1);
 

for (i in 1:nrow(bfmat)){
#for (i in 1:5){	
	for (j in 1:nrow(bfmat)){
		if (!is.na(bfmat[i,j])){
			xval <- as.numeric(rownames(bfmat)[i]);
			yval <- as.numeric(colnames(bfmat)[j])
			xco <- c(xval, xval, xval+1, xval+1);
			yco <- c(yval, yval+1, yval+1, yval);
			polygon(x=xco, y=yco, lwd=0.8, col=cmat[i, j], border=cmat[i,j])
		}
	}
}

axis(1, at=seq(-10, 100, by=10), cex.axis=1.2);
axis(2, at=seq(-10, 100, by=10), las=1, cex.axis=1.2);
mtext(side=1, text="Macroevolutionary regimes, numerator", cex=1.4, line=3.5);
mtext(side=2, text="Macroevolutionary regimes, denominator", cex=1.4, line=3.5);

# Add the color bar
cb <- getColorBar();

plot.new();
par(mar=c(7,1,1,7));

plot.window(xlim=c(0, 2.5), ylim=c(-8,8));

for (i in 2:nrow(cb)){
	xv <- c(0,0,1,1);
	yv <- c(cb$cuts[c(i-1,i,i,i-1)]);
	polygon(xv, yv, col=cb$cols[i]);
	
}

axis(4, at=seq(-5, 5, by=2), las=1, cex.axis=2,pos=1.5);
mtext(side=4, at=0, text="(Log) Bayes factor", cex=1.5)

 


