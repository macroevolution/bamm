
library(BAMMtools);


#
# We will use the pre-loaded event data for the whales:
#

data(events.whales); # here is the event data
data(whales); # here is the whale tree

# Process the event data in a 'bammdata' object:
bammdata <- getEventData(whales, events.whales, burnin=0.1);

class(bammdata);  # check: should be 'bammdata'

# We will plot the 1st, 30th, and 40th samples from the posterior.

ixx <- rep(c(10, 30, 40), 3);

# Set up the plot windows
par(mar=numeric(4));
quartz.options(height=10, width=10, dpi=72); #This line is OSX specific
plot.new();
par(mfrow=c(3,3));
 
# Here we set up the color schemes 
colschemes <- list();
colschemes[1:3] <- 'temperature';
colschemes[4:6] <- 'Spectral';
colschemes[7:9] <- list(c('blue', 'gray', 'red'));

for (i in 1:length(ixx)) {
#for (i in 1:3){
	
	index <- ixx[i];
	
	eventsub <- subsetEventData(bammdata, index=index);
	plot.bammdata(eventsub, method='polar', pal= colschemes[[i]], par.reset=FALSE);
	addBAMMshifts(eventsub, method='polar', index=1, col='white', bg='black', cex=4, par.reset=FALSE);
}
