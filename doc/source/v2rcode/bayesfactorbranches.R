
# 
data(whales, prior.whales, events.whales)
edata <- getEventData(whales, events.whales, burnin=0.1)

plot(whales, show.tip.label=F)
margprobs <- marginalShiftProbsTree(edata)

prior <- getBranchShiftPriors(whales, prior.whales)
bfmat <- bayesFactorBranches(edata, prior)

index <- which(margprobs$edge.length > 0.05)
margprobs$edge[index,2]
prior$edge.length[index]
 
 



x1 <- c(22.71, 16.13)
y1 <- c(56.96, 66.26)
 
 
# This part is just to locate...
ew <- rep(1, length(whales$edge.length))
ew[index] <- 4
ec <- rep('black', length(whales$edge.length))
ec[index] <- 'red'

plot(whales, show.tip.label=F, no.margin=T, edge.width=ew, edge.color=ec, x.lim=c(0, 45))

ww <- locator() 
# first pick the text, then the line:
arrows(x0 = ww$x[2], x1 = ww$x[3] , y0 = ww$y[2], y1 =ww$y[3], length=0.1)
text(x=ww$x[1], y=ww$y[1],label='0.37', cex=1.5, col='red', font=4) 

zz <- locator() 
# first pick the text, then the line:
arrows(x0 = zz$x[2], x1 = zz$x[3] , y0 = zz$y[2], y1 = zz$y[3], length=0.1)
text(x=zz$x[1], y=zz$y[1],label='0.55', cex=1.5, col='red', font=4) 

yy <- locator() 
# first pick the text, then the line:
arrows(x0 = yy$x[2], x1 = yy$x[3] , y0 = yy$y[2], y1 = yy$y[3], length=0.1)
text(x=yy$x[1], y=yy$y[1],label='0.06', cex=1.5, col='red', font=4) 


####### Now for png plots:::
png(file = 'bayesfactorbranches1.png')
ew <- rep(1, length(whales$edge.length))
ew[index] <- 4
ec <- rep('black', length(whales$edge.length))
ec[index] <- 'red'

plot(whales, show.tip.label=F, no.margin=T, edge.width=ew, edge.color=ec, x.lim=c(0, 45))
 
# first pick the text, then the line:
arrows(x0 = ww$x[2], x1 = ww$x[3] , y0 = ww$y[2], y1 =ww$y[3], length=0.1)
text(x=ww$x[1], y=ww$y[1],label='0.37', cex=1.5, col='red', font=4) 

 # first pick the text, then the line:
arrows(x0 = zz$x[2], x1 = zz$x[3] , y0 = zz$y[2], y1 = zz$y[3], length=0.1)
text(x=zz$x[1], y=zz$y[1],label='0.55', cex=1.5, col='red', font=4) 

 # first pick the text, then the line:
arrows(x0 = yy$x[2], x1 = yy$x[3] , y0 = yy$y[2], y1 = yy$y[3], length=0.1)
text(x=yy$x[1], y=yy$y[1],label='0.06', cex=1.5, col='red', font=4) 

dev.off()


png(file = 'bayesfactorbranches1.png')
ew <- rep(1, length(whales$edge.length))
ew[index] <- 4
ec <- rep('black', length(whales$edge.length))
ec[index] <- 'red'

plot(whales, show.tip.label=F, no.margin=T, edge.width=ew, edge.color=ec, x.lim=c(0, 45))
 
# first pick the text, then the line:
arrows(x0 = ww$x[2], x1 = ww$x[3] , y0 = ww$y[2], y1 =ww$y[3], length=0.1, lwd=2)
text(x=ww$x[1], y=ww$y[1],label='0.37', cex=1.5, col='red', font=4) 

 # first pick the text, then the line:
arrows(x0 = zz$x[2], x1 = zz$x[3] , y0 = zz$y[2], y1 = zz$y[3], length=0.1, lwd=2)
text(x=zz$x[1], y=zz$y[1],label='0.55', cex=1.5, col='red', font=4) 

 # first pick the text, then the line:
arrows(x0 = yy$x[2], x1 = yy$x[3] , y0 = yy$y[2], y1 = yy$y[3], length=0.1, lwd=2)
text(x=yy$x[1], y=yy$y[1],label='0.06', cex=1.5, col='red', font=4) 

dev.off()


png(file = 'bayesfactorbranches2.png')
ew <- rep(1, length(whales$edge.length))
ew[index] <- 4
ec <- rep('black', length(whales$edge.length))
ec[index] <- 'darkgreen'

plot(whales, show.tip.label=F, no.margin=T, edge.width=ew, edge.color=ec, x.lim=c(0, 45))
 
# first pick the text, then the line:
arrows(x0 = ww$x[2], x1 = ww$x[3] , y0 = ww$y[2], y1 =ww$y[3], length=0.1, lwd=2)
text(x=ww$x[1], y=ww$y[1],label='0.008', cex=1.5, col='darkgreen', font=4) 

 # first pick the text, then the line:
arrows(x0 = zz$x[2], x1 = zz$x[3] , y0 = zz$y[2], y1 = zz$y[3], length=0.1, lwd=2)
text(x=zz$x[1], y=zz$y[1],label='0.001', cex=1.5, col='darkgreen', font=4) 

 # first pick the text, then the line:
arrows(x0 = yy$x[2], x1 = yy$x[3] , y0 = yy$y[2], y1 = yy$y[3], length=0.1, lwd=2)
text(x=yy$x[1], y=yy$y[1],label='0.025', cex=1.5, col='darkgreen', font=4) 

dev.off()


png(file = 'bayesfactorbranches3.png')
ew <- rep(1, length(whales$edge.length))
ew[index] <- 4
ec <- rep('black', length(whales$edge.length))
ec[index] <- 'blue'

plot(whales, show.tip.label=F, no.margin=T, edge.width=ew, edge.color=ec, x.lim=c(0, 45))
 
# first pick the text, then the line:
arrows(x0 = ww$x[2], x1 = ww$x[3] , y0 = ww$y[2], y1 =ww$y[3], length=0.1, lwd=2)
text(x=ww$x[1], y=ww$y[1],label='69', cex=1.5, col='blue', font=4) 

 # first pick the text, then the line:
arrows(x0 = zz$x[2], x1 = zz$x[3] , y0 = zz$y[2], y1 = zz$y[3], length=0.1, lwd=2)
text(x=zz$x[1], y=zz$y[1],label='819', cex=1.5, col='blue', font=4) 

 # first pick the text, then the line:
arrows(x0 = yy$x[2], x1 = yy$x[3] , y0 = yy$y[2], y1 = yy$y[3], length=0.1, lwd=2)
text(x=yy$x[1], y=yy$y[1],label='2.6', cex=1.5, col='blue', font=4) 

dev.off()

### 
png(file = 'bayesfactorbranches4.png')
ew <- rep(1, length(whales$edge.length))
ew[index] <- 4
ec <- rep('black', length(whales$edge.length))
ec[index] <- 'red'
plot.phylo(bf, show.tip.label=F, edge.color=ec, edge.width=ew)
lines(x=c(200, 300), y=c(20,20), lwd=3, col='blue')
text(x=300, y=25, label='BF = 100', cex=1.5, font=3)
dev.off()























