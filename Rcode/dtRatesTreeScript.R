library(ape); library(gplots);
source('plot.dtrates.R');
#In terminal you will have to run R CMD SHLIB dtrates.c
#Hopefully you will have the right compiler in your path; if so,
#This will create two files, and dtrates.so is the one you want
dyn.load('dtrates.so');

phy = read.tree('whaletree.tre');
fn = 'beventdata.txt';

x = getEventData(phy, fn, burnin = 0.5);
x = dtRates(x, 0.01);

p = plot.dtrates(x,method='phylogram',lwd=3,palette='temperature',ncolors=64);
p = plot.dtrates(x,method='polar',lwd=3,vtheta=5,rbf=0.001,palette='temperature',ncolors=64);

p = plot.dtrates(x,method='phylogram',lwd=3,palette='diverging',ncolors=64);
p = plot.dtrates(x,method='polar',lwd=3,vtheta=5,rbf=0.001,palette='diverging',ncolors=64);