library(ape); library(gplots);
source('dtRatesTreeFxns.R');

phy = read.tree('whaletree.tre');
fn = 'beventdata.txt';

x = getEventData(phy, fn, burnin = 0.5);
x = dtRates(x, 0.01);

polartree(x,5,lwd=3,colorize=TRUE);