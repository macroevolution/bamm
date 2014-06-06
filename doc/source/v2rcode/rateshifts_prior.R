library(BAMMtools)
data(prior.whales)

png(file = 'rateshifts_prior.png', width=450, height=450);

#quartz.options(height=6, width=6)
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0,18), ylim=c(0,0.45))

lines(prior.whales$N_shifts, prior.whales$prob, lwd=1.5, col='gray30')
points(prior.whales$N_shifts, prior.whales$prob, pch=21, bg='red', cex=1.5)
axis(1, at=seq(-2, 18, by=2), cex.axis=1.2)
axis(2, at=seq(-.2, 0.6, by=0.2), cex.axis=1.2, las=1)
mtext("Probability", side=2, line=3.5, cex=1.5)

mtext("Number of shifts", side=1, line=3, cex=1.5)

dev.off();

##### 

xx <- read.csv('prior_probs_0.1.txt')

png(file = 'rateshifts_prior_0.1.png', width=450, height=450);

#quartz.options(height=6, width=6)
plot.new()
par(mar=c(5,5,1,1))
plot.window(xlim=c(0,18), ylim=c(0,0.45))

lines(xx$N_shifts, xx$prob, lwd=1.5, col='gray30')
points(xx$N_shifts, xx$prob, pch=21, bg='red', cex=1.5)
axis(1, at=seq(-2, 18, by=2), cex.axis=1.2)
axis(2, at=seq(-.2, 0.6, by=0.2), cex.axis=1.2, las=1)
mtext("Probability", side=2, line=3.5, cex=1.5)

mtext("Number of shifts", side=1, line=3, cex=1.5)

dev.off();



