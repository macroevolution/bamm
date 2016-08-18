

require(deSolve)
require(Rcpp)

sourceCpp("MonteCarloSimulator.cpp")

# functions required for the analyses:
source("MEA_analysis_functions.R")

pars1 <- c(0.001, 0.0009, 0.001, 50, 50)
names(pars1) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")

pars0 <- c(1, 1.1, 0.001, 50, 50)
names(pars0) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")

# initial segment of x0 = 1, rate shift, then 49 time units to present
#    with pars1
tvec <- c(0, 1, 50)
 

compute_MEA_probability(pars0, pars1, tvec)

compute_BAMM_probability(pars0 , pars1, tvec)


############################################################

sourceCpp("MonteCarloSimulator.cpp") 
# this file from MEA Dryad, supplementary_data/code/monte_carlo_simulation/

# add Moore et al likelihood functions
source("likelihoodFunctions.R") 
#this file from MEA Dryad, supplementary_data/code/monte_carlo_simulation/

# our file, contains wrapper for MEA functions
source("MEA_analysis_functions.R")

MEA_prob_PNAS_implementation(pars0, pars1, tmax=50, tshift = 1) 


############################################################


# looping over some values

rel_extinction <- 0.999

lambda_vec <- seq(0.001, 3, by=0.1)
mu_vec <- lambda_vec * rel_extinction

mea_probs <- numeric(length(lambda_vec))
bamm_probs <- numeric(length(mu_vec))

# keeping other params same as previous
pars1 <- c(0.001, 0.0009, 0.001, 50, 50)
names(pars1) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")
tvec <- c(0, 1, 50)

for (i in 1:length(lambda_vec)){
	cat(i, "\n")
	ipars <- c(lambda = lambda_vec[i], mu=mu_vec[i], pars1[2:5])
	
	mea_probs[i] <- compute_MEA_probability(ipars, pars1, tvec)
	bamm_probs[i] <- compute_BAMM_probability(ipars, pars1, tvec)
}


########### plotting the exercise above ##################

quartz.options(height=6, width=6)
plot.new()
par(oma=c(1,1,1,1))
 

plot.new()
par(mar=c(6,6,1,1))
plot.window(xlim=c(-0.1, 3.1), ylim=c(0,8))
lines(x=c(-.5, 3), y=c(1,1), lwd=3, col="gray60")
points(lambda_vec, mea_probs, pch=21, bg="red", cex=1.3)
points(lambda_vec, bamm_probs, pch=21, bg="blue", cex=1.3)

axis(2, at=seq(-1, 8, by=1), las=1)
axis(1, at=seq(-0.5, 3, by=0.5), las=1)

mtext(side=1, text = "Speciation rate", cex=1.5, line=3)
mtext(side=2, text = "Probability", cex=1.5, line=3)

###########################################
 

mu_vec <- seq(0.15, 2, by=0.1)

mea_probs_highE <- numeric(length(mu_vec))
bamm_probs_highE <- numeric(length(mu_vec))

# keeping other params same as previous
pars1 <- c(0.001, 0.0009, 0.001, 50, 50)
names(pars1) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")
tvec <- c(0, 1, 50)

for (i in 1:length(mu_vec)){
	cat(i, "\n")
	ipars <- c(lambda = 0.1, mu=mu_vec[i], pars1[2:5])
	
	mea_probs_highE[i] <- compute_MEA_probability(ipars, pars1, tvec)
	bamm_probs_highE[i] <- compute_BAMM_probability(ipars, pars1, tvec)
}

########### plotting the exercise above for extinction ##################

quartz.options(height=6, width=10)
plot.new()
par(oma=c(1,1,1,1))
par(mfcol=c(1,2))

plot.new()
par(mar=c(6,6,1,1))
plot.window(xlim=c(-0.1, 2.1), ylim=c(0, 400))
lines(x=c(-.5, 3), y=c(1,1), lwd=3, col="gray60")
points(mu_vec, mea_probs_highE, pch=21, bg="red", cex=1.3)
points(mu_vec, bamm_probs_highE, pch=21, bg="blue", cex=1.3)

axis(2, at=seq(-50, 400, by=50), las=1)
axis(1, at=seq(-0.5, 3, by=0.5), las=1)

mtext(side=1, text = "Extinction rate", cex=1.5, line=3)
mtext(side=2, text = "Probability", cex=1.5, line=3)

# on a log-scale:

plot.new()
par(mar=c(6,6,1,1))
plot.window(xlim=c(-0.1, 2.1), ylim=c(-1, 3.1))
lines(x=c(-.5, 3), y=c(0,0), lwd=3, col="gray60")
points(mu_vec, log10(mea_probs_highE), pch=21, bg="red", cex=1.3)
points(mu_vec, log10(bamm_probs_highE), pch=21, bg="blue", cex=1.3)

axis(2, at=seq(-2, 3, by=1), las=1, labels=c(0.01, 0.1, 1, 10, 100, 1000))
axis(1, at=seq(-0.5, 3, by=0.5), las=1)

mtext(side=1, text = "Extinction rate", cex=1.5, line=3)
mtext(side=2, text = "Probability", cex=1.5, line=3)



########### alternative parameterization ##################
 
pars1 <- c(0.001, 0.0009, 0.001, 50, 50)
names(pars1) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")

pars0 <- c(0.1, 0.08, 0.001, 50, 50)
names(pars0) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")
 
tvec <- c(0, 1, 50)
 
x <- replicate(1e5, SimulateCPBDP(50, pars0["lambda"], pars0["mu"], pars0["shiftrate"], pars0["prior_lambda"], pars0["prior_mu"])) 
mean(x < 50)
# = 0.715, far from potential boundary effects at 0 or 1
 
compute_MEA_probability(pars0, pars1, tvec)

compute_BAMM_probability(pars0 , pars1, tvec)

 
MEA_prob_PNAS_implementation(pars0, pars1, tmax=50, tshift = 1) 


########### Illustration using analytical equations ##############################
########### This model disallows rate shifts on extinct branches ##################



pars1 <- c(0.001, 0.0009, 0.001, 50, 50)
names(pars1) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")

pars0 <- c(0.1, 0.08, 0.001, 50, 50)
names(pars0) <- c("lambda", "mu", "shiftrate", "prior_lambda", "prior_mu")
 
tvec <- c(0, 1, 50)

compute_MEA_analytical_noExtinctShifts(pars0, pars1, tvec)


