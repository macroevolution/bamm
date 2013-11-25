library(ape); library(gplots);
source('plot.dtrates.R');

phy = read.tree('whaletree.tre');
fn = 'beventdata.txt';

x = getEventData(phy, fn, burnin = 0.5);
x = dtRates(x, 0.01);

plot.dtrates(x,method='phylogram',lwd=3,palette='temperature',ncolors=64);
plot.dtrates(x,method='polar',lwd=3,vtheta=5,rbf=0.001,palette='temperature',ncolors=64);