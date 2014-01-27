
library(BAMMtools);
data(primates);
data(events.primates);

ed <- getEventData(primates, events.primates, burnin=0.2, type='trait');

makeTransparent<-function(someColor, alpha=100)
{
  	newColor<-col2rgb(someColor)
  	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

graynew <- makeTransparent('gray40', alpha=75);

quartz.options(height=6, width=10, dpi=72); # on mac osx 
par(oma=c(0,0,0,0));
par(mar=c(1,1,1,1));
par(mfrow=c(1,2));

plot.bammdata(ed, lwd=2.5, method='polar', pal = c("darkgreen", "yellow2", "red"));

plot.bammdata(ed, lwd=2.5, method='polar', pal = c("darkgreen", "yellow2", "red"));

for (i in 1:400){
	addBAMMshifts(ed, method='polar', index=i, cex=1.5, pch=19, col =graynew);	
}




 