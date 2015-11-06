# This supplement contains functions and scripts 
# to repeat simple calculations detailed on the website
#
# The initial block of functions are required for the simulation script
# and include implementations of the D(t) and E(t) equations for the constant
# rate birth-death process.
# 
# functions "simulateHistory_typeI" etc simulate a specified diversification process
# and return TRUE if the process is extinct in the present, and FALSE if it 
# does not go extinct.
#
# E_func and D_func are the D(t) and E(t) equations from Rabosky 2014 (PLoS ONE)
 

################################################
################################################
############ Functions

 
# likelihood functions for the constant-rate birth-death process
#    used to compute speciation and extinction probabilities
#    on individual branches.

E_func <- function(lam, mu, E0, dt) {
	
	num <- (1 - E0) * (lam - mu);
	denom <- (1 - E0) * lam - (mu - lam * E0) * exp(-(lam - mu) * dt);
	
	return( 1 - num / denom)	
 
}

D_func <- function(lam, mu, E0, D0, dt) {
	
	r <- lam - mu
	num <- (D0 * r^2) * exp((-r) * dt);
	denom <- ( (lam - lam * E0 + exp((-r) * dt) * (lam * E0 - mu)  )  ) ^ 2;
	return(num / denom) 
 
}


betaFx <- function(lambda, mu, dt){
	rr <- lambda - mu
	eps <- mu / lambda
	beta <- ( (exp(rr * dt) - 1 ) / (exp(rr * dt) - eps))	
	return(beta)
}


getTotalProgeny <- function(n, lambda, mu, dt){
	beta <- betaFx(lambda, mu, dt)
	xvec <- numeric(n)
	
	if (n == 0){
		return(0)
	}
	
	# the extinction probability, after Raup 1985
	p0 <- (mu / lambda) * beta
	for (i in 1:length(xvec)) {
		isExtinct <- runif(1) <= p0
		if (isExtinct){
			xvec[i] <- 0
		}else{
			# if not extinct, draw from shifted geometric distribution
			xvec[i] <- rgeom(1, prob=(1 - beta)) + 1
		}
	}
	return(sum(xvec))
}

 
# Type 1. single lineage undergoes the shift
# The process only survives if a lineage 
# survives to the present having undergone all of the mapped
# diversification histories

simulateHistory_type1 <- function(x){
 
	dt <- x$t_end - x$t_start
 
    # each interval just starts with 1 lineage
    # this is all that matters, since only a single 
    # lineage is chosen for the rate shift
    
	for (i in 1:nrow(x)){
		NN <- getTotalProgeny(1, x$lambda[i], x$mu[i], dt[i])
		if (NN == 0){
			return(TRUE)
		}
	}
    return(FALSE)
}

# Type 2. This function simulates an instance of a process and returns TRUE 
# if the process goes extinct, false otherwise.
# Model assumes a discrete diversification history
# along branches, with instantaneous changepoints, such that all 
# lineages shift instantly to a new diversification rate 
# at each shift point.
# This is the BAMM process.

simulateHistory_type2 <- function(x){
 
	dt <- x$t_end - x$t_start

	NN <- 1
    	
	for (i in 1:nrow(x)){
		NN <- getTotalProgeny(NN, x$lambda[i], x$mu[i], dt[i])
		if (NN == 0){
			return(TRUE)
		}
	}
    return(FALSE)
}

# Type 3. this model (type III) simulates an instance of a process 
# but where extinction only happens when all lineages generated
# go extinct. Rate shifts are assumed to happen to just a single lineage
# but the process can still survive as a whole even if the surviving lineage
# does not include the fully mapped diversification history 

simulateHistory_type3 <- function(x){
 
	dt <- x$t_end - x$t_start
 
	mm <- matrix(0, nrow=nrow(x), ncol=nrow(x))
    
  	 	# generate initial progeny vector
	mm[1,1] <- getTotalProgeny(1, x$lambda[1], x$mu[1], dt[1])
	start_vec <- numeric(nrow(mm))  	
    	
    	
	for (i in 2:ncol(mm)){
    	# initialize:
    	
		if (mm[i-1,i-1] == 0){
			break;
    	}else{
			start_vec <- mm[,i-1]
			
			# choose an individual to the the next "type"
			chosen <- sample(1:length(start_vec), size=1, prob = start_vec / sum(start_vec))
			start_vec[chosen] <- start_vec[chosen] - 1
			start_vec[i] <- 1
 
		}

		for (j in 1:nrow(mm)){
			if (start_vec[j] > 0){
				mm[j,i] <- getTotalProgeny(start_vec[j], x$lambda[j], x$mu[j], dt[j])
				if (i == ncol(mm) & mm[j,i] > 0){
					# process has not gone extinct by end of simulation
					return(FALSE)
				} 				
			}
  		
				
		}
	}
 	return(TRUE)
}


# Compute extinction probability along a single branch
# using the computations as currently implemented in BAMM

computeExtinctionBAMM <- function(x){
	
	evec <- numeric(nrow(x))
	for (i in nrow(x):1){
		#cat(x$t_start[i], "\n")
		dt <- x$t_end[i] - x$t_start[i]
		
		E0 <- 0
		if (i != nrow(x)){
			E0 <- evec[i+1] 
		}
		
		evec[i] <- E_func(x$lambda[i], x$mu[i], E0, dt)
	
	}
	return(evec[1])
}


# Compute extinction probability along a single branch
# using the "recomputing" algorithm.

computeExtinctionRecomputed <- function(x){
	
	evec <- numeric(nrow(x))
	
	maxT <- max(x$t_end) 
	
	for (i in nrow(x):1){
		#cat(x$t_start[i], "\n")
		dt <- x$t_end[i] - x$t_start[i]
		
		E0 <- 0
		if (i != nrow(x)){
			E0 <- E_func(x$lambda[i], x$mu[i], 0, (maxT - x$t_end[i]))
		}
		
		evec[i] <- E_func(x$lambda[i], x$mu[i], E0, dt)
	
	}
	return(evec[1])	
	
}

# 
simulateExtinctionProb <- function(x, nsims=2000, type = 1){
	
	res <- numeric(nsims)
	for (i in 1:nsims){
		if (type == 1){
			res[i] <- simulateHistory_type1(x)
		}else if (type == 2){
			res[i] <- simulateHistory_type2(x)
		}else if (type == 3){
			res[i] <- simulateHistory_type3(x)
		}else{
			stop("invalid option")
		}	
	}
	return(sum(res == TRUE) / length(res))		
	
}



######################
######################
######################

# HEre we simulate extinction probabilities and compare them between BAMM
#      and recomputed implementations


NSIMS <- 500
ESIMS <- 5000

x <- numeric(NSIMS)
res <- data.frame(sim1=x, sim2 = x, sim3=x, bamm=x, recomputed=x, lam1=x, lam2=x, mu1=x, mu2=x)

for (i in 1:NSIMS){
	cat(i, '\n')
	eps <- runif(2, 0, 1)
	lambda <- runif(2, 0, 0.05)
	mu <- lambda * eps
 
	
	t_start <- c(0, 20)
	t_end <- c(20, 100)	
	dff <- data.frame(t_start, t_end, lambda, mu)
	
	tmp <- numeric(ESIMS)
 
	res$sim1[i] <- simulateExtinctionProb(dff, type = 1, nsims = ESIMS)
	res$sim2[i] <- simulateExtinctionProb(dff, type = 2, nsims = ESIMS)
	res$sim3[i] <- simulateExtinctionProb(dff, type = 3, nsims = ESIMS)
	res$bamm[i] <-  computeExtinctionBAMM(dff)
	res$recomputed[i] <- computeExtinctionRecomputed(dff)

	res$lam1[i] <- dff$lambda[1]
	res$lam2[i] <- dff$lambda[2]
	res$mu1[i] <- dff$mu[1]
	res$mu2[i] <- dff$mu[2]
 
}

write.table(res, file = "extinction_simulation_results.csv", sep=",", quote=F, row.names=F)


 ############# plotting results
 
xx <- read.csv("extinction_simulation_results.csv", header=T) 
 
plotSetup <- function(){
	plot.new()
	par(mar=c(6,6,1,1))
	plot.window(xlim=c(0,1), ylim=c(0,1))
	lines(x=c(0,1), y=c(0,1), lwd=2, col="gray50")
	axis(1, at=seq(-0.2, 1, by=0.2), cex.axis=1.5)
	axis(2, at=seq(-0.2, 1, by=0.2), las=1, cex.axis=1.5)
}
 
 
 
#quartz.options(height = 4.5, width = 10)
png(height = 450, width=1000, file = "x_extinctionprobs_bamm.png")
plot.new()
par(mfrow=c(1,3))
par(oma=c(1,1,3,1))

plotSetup()
points(xx$sim1, xx$bamm, pch=19, col="blue", cex=0.8)
mtext(side = 2, text = "BAMM extinction probability", line = 4, cex = 2 ,font=2)
mtext(side = 3, text = "Type 1 shift", cex = 2.2, font=2)

plotSetup()
points(xx$sim2, xx$bamm, pch=19, col="blue", cex=0.8)
mtext(side = 1, text = "Simulated extinction probability under focal process", line = 4.5, cex = 2, font=2)
mtext(side = 3, text = "Type 2 shift", cex = 2.2, font=2)

plotSetup()
points(xx$sim3, xx$bamm, pch=19, col="blue", cex=0.8)
mtext(side = 3, text = "Type 3 shift", cex = 2.2, font=2)
dev.off()


 #quartz.options(height = 4.5, width = 10)
png(height = 450, width=1000, file = "x_extinctionprobs_recomputed.png")
plot.new()
par(mfrow=c(1,3))
par(oma=c(1,1,3,1))

plotSetup()
points(xx$sim1, xx$recomputed, pch=19, col="red", cex=0.8)
mtext(side = 2, text = "RECOMPUTED extinction probability", line = 4, cex = 2 ,font=2)
mtext(side = 3, text = "Type 1 shift", cex = 2.2, font=2)

plotSetup()
points(xx$sim2, xx$recomputed, pch=19, col="red", cex=0.8)
mtext(side = 1, text = "Simulated extinction probability under focal process", line = 4.5, cex = 2, font=2)
mtext(side = 3, text = "Type 2 shift", cex = 2.2, font=2)

plotSetup()
points(xx$sim3, xx$recomputed, pch=19, col="red", cex=0.8)
mtext(side = 3, text = "Type 3 shift", cex = 2.2, font=2)
dev.off()


 











