
 
# likelihood functions for the constant-rate birth-death process
#    used to compute speciation and extinction probabilities
#    on individual branches.
# These functions provide the analytical solutions to D(t) and E(t)
# for arbitrary initial speciation and extinction probabilities (D0, E0)
#
#
#  
# E_func(...): Analytical solution to E(t)
#          args: lam = speciation rate
#                mu  = extinction rate
#                E0  = initial extinction probability
#                dt  = time interval
# 
# return value: extinction prob after dt 
 
E_func <- function(lam, mu, E0, dt) {
	
	num <- (1 - E0) * (lam - mu);
	denom <- (1 - E0) * lam - (mu - lam * E0) * exp(-(lam - mu) * dt);
	
	return(as.vector( 1 - num / denom))
 
}

#
# D_func(...): Analytical solution to D(t)
#          args: lam = speciation rate
#                mu  = extinction rate
#                D0  = initial prob of data
#                E0  = initial extinction probability
#                dt  = time interval
# 
# return value: extinction prob after dt 
D_func <- function(lam, mu, E0, D0, dt) {
	
	r <- lam - mu
	num <- (D0 * r^2) * exp((-r) * dt);
	denom <- ( (lam - lam * E0 + exp((-r) * dt) * (lam * E0 - mu)  )  ) ^ 2;
	return(as.vector(num / denom)) 
 
}
 
 
 
 
speciationFunctionFast <- function(t, y, parms, evec, e_init = 0){
	with (as.list(c(y, parms)), {  		
		
		E <- e_init 
	 	if (t > 0){
	 	 	E <- sum(evec < t) / length(evec)
	
	 	}
 
 		# process should include shiftrate in first term -(lambda + mu + shiftrate)
 		# MEA do not, but has relatively small effect here and does 
 		# not cause serious problems 
		#dD <- - (lambda  + mu  + shiftrate) * D_t + 2*lambda * E * D_t
		
		# equation exactly as used by MEA
		dD <- - (lambda  + mu  ) * D_t + 2*lambda * E * D_t
		list(dD)
	})
} 

compute_MEA_probability <- function(pars0, pars1, tvec){
	
	# pars0 is initial process at root
	# pars1 is derived process, after rate shift
	
	dt1 <- tvec[3] - tvec[2]
	dt0 <- tvec[2] - tvec[1]
	
	# note: dt0 is also the shift time since 
	#       we have reset vector to start at t_init = 0
	
	total_age <- tvec[3] - tvec[1]
 
 	# generate extinction probability vectors
 	exprob1 <- replicate(100000, SimulateCPBDP(dt1, pars1['lambda'], pars1['mu'], 
 								pars1['shiftrate'], pars1['lambda_prior'], pars1['mu_prior']))
 	exprob0 <- replicate(100000, SimulateCPBDP(total_age, pars0['lambda'], pars0['mu'], 
 								pars0['shiftrate'], pars0['lambda_prior'], pars0['mu_prior']))
		
	
	# initial probability of data (= 1)
	cur_state <- c(D_t = 1)
	
	d1 <- ode(y = cur_state, times = c(0,dt1), func = speciationFunctionFast, parms = pars1, evec = exprob1)
	

	
 	# this value of d1 becomes input for second interval
 	# as does exprob2 computed between shifttime and present
 	cur_state <- c(D_t = as.numeric(d1[2,2]))
	
	# intial prob of extinction for the integration, given initial process with pars0
	e_init_seg0 <- sum(exprob0 < dt1) / length(exprob0)
	
	d2 <- ode(y = cur_state, times = c(0, dt0), func = speciationFunctionFast, parms = pars0, evec = exprob0, e_init = e_init_seg0)
	
	exprob_at_root <- sum(exprob0 < total_age) / length(exprob0)
	
	# probability of data
	#  conditioned on single lineage survival:
	prob <- as.numeric(d2[2,2]) / (1 - exprob_at_root)
	return(prob)
}


compute_BAMM_probability <- function(pars0, pars1, tvec){
	
	dt1 <- tvec[3] - tvec[2]
	dt0 <- tvec[2] - tvec[1]
	
	# note: dt0 is also the shift time since 
	#       we have reset vector to start at t_init = 0
	
	total_age <- tvec[3] - tvec[1]
	
	d1 <- D_func(pars1["lambda"], pars1["mu"], E0 = 0, D0 = 1, dt = dt1)
	e1 <- E_func(pars1["lambda"], pars1["mu"], E0 = 0, dt = dt1) 
	d2 <- D_func(pars0["lambda"], pars0["mu"], E0 = as.numeric(e1), D0 = as.numeric(d1), dt= dt0)
	e2 <- E_func(pars0["lambda"], pars0["mu"], E0 = as.numeric(e1), dt = dt0)
	
	bamm_prob <- as.numeric(d2 / (1 - e2))
	
	return(bamm_prob)
	
}

 
# These functions are wrappers which use the code from Moore et al to compute 
# the probability of a single branch with an augmented history (eg, mapped rate shift)
 ##########################################################################
# calcBranchProb() is a simple wrapper to reformat parameter combinations and 
# use the Moore et al. computeBranchProbability() function to calculate the 
# probability of a single *terminal* branch with a shift
	
calcBranchProb <- function(mus=c(0.25,0.25), # root and derived extinction rates
							lambdas=c(1,1),  # root and derived speciation rates
							ShiftT=0.1, 	 # time from root for shift
							endT=1, 		 # total length of the edge
							grainLen=1e4,	 # grain of the intervals 
							nreps=1e4, 		 # number of simulations for extinction
							transition_rate=0.1, # transition rate to new processes
							lambda_prior=1,  # prior on speciation for new processes
							mu_prior=1)	{  # prior on extinction for new processes
								
	# calcBranchProb() reformats arbitrary data into a data.frame object that 
	# computeBranchProbability() can use as an argument
	
	# simulated_times_root is a vector of extinction times for the root process, 
	# to determine the probability it will go extinct at an arbitrary time
	simulated_times_root <- replicate(nreps, SimulateCPBDP(endT, lambdas[1], mus[1], 
												transition_rate, lambda_prior, mu_prior))
												
	# e_prob_at_root is the probability of the root regime going extinct before endT											
	e_prob_at_root <- mean(simulated_times_root < endT)
	
	# simulated_times_shift is a vector of extinction times for the shifted process
	simulated_times_shift <- replicate(nreps, SimulateCPBDP(endT-ShiftT, lambdas[2], mus[2], 
												transition_rate, lambda_prior, mu_prior)) 	
	
	# There may be an implementation issue in the original MEA code 
	# with these_intervals without adding ShiftT; we use the original in this implementation to 
	# avoid modifications to MEA code, but the snippet below can be uncommented to test this
	#
	# The shifted process is endT-ShiftT in length, but it occurs on intervals from
	# ShiftT *to* endT. Hence the need to add ShiftT to yield extinction times *from root*
	#simulated_times_shift <- replicate(nreps, SimulateCPBDP(endT-ShiftT, lambdas[2], mus[2], 
	#											transition_rate, lambda_prior, mu_prior)) + ShiftT
	#
	# 


	# computeBranchProbability() requires a dataframe with the specation and extinction rates,
	# as well as the simulated times of extinction. Columns unused by computeBranchProbability()
	# are filled with NA
	processes <- data.frame(
					event_time=c(NA, NA), 
					node=c(NA, NA), 
					node_time=c(NA, NA), 
					speciation=lambdas, 
					extinction=mus, 
					simulated_times=I(list(simulated_times_root, simulated_times_shift)))

	# these_intervals are the delta-T levels computeBranchProbability calculates values over
	these_intervals <- seq(0, endT, length.out=grainLen)
	
	# these_processes is a vector of process-to-interval mappings, so computeBranchProbability()
	# know what process governs each small segment of the branch
	these_processes <- rep(1, length(these_intervals))
	these_processes[which(these_intervals > ShiftT)] <- 2
	branch_intervals <- list(these_intervals=these_intervals, these_processes=these_processes)

	# The above code reformats parameters into an object that can be passed to the
	# Moore et al. (2016) function from likelihood.R. This outputs 
	# the probability of an individual branch
	out <- computeBranchProbability(1, branch_intervals, processes)
	return(c(out, e_prob_at_root))
}

# MEA_prob_PNAS_implementation() is another wrapper to take parameter vectors and run the calcBranchProb() function
# pars0 is a named vector of lambda, mu, shiftrate, prior_lambda & prior_mu for the root regime
# pars1 is the same as pars0 but for the derived, shifted regime
MEA_prob_PNAS_implementation <- function(pars0, pars1, tmax, tshift, Nreps=1e5)	{
	# This calculates the probability of a branch with a mapped shift
	Prob <- calcBranchProb(mus=c(pars0["mu"], pars1["mu"]), 
			lambdas=c(pars0["lambda"], pars1["lambda"]), 
			ShiftT=tshift, endT=tmax, 
			transition_rate=pars0["shiftrate"], 
			lambda_prior=pars0["prior_lambda"], 
			mu_prior=pars0["prior_mu"], nreps=Nreps)
	
	# out is the probability of the branch (Prob[1]) divied by the probability the root process
	# survived to the present (1 - Prob[2])	
	out <- Prob[1] / (1 - Prob[2])
	return(out)
}

 

# this function computes MEA-type probability ("recompute") but
# disallows shifts on extinct / unobserved branches

compute_MEA_analytical_noExtinctShifts <- function(pars0, pars1, tvec){
	
	dt1 <- tvec[3] - tvec[2]
	dt0 <- tvec[2] - tvec[1]
	total_age <- tvec[3] - tvec[1]	
	
	# data probability for initial segment
	d1 <- D_func(pars1["lambda"], pars1["mu"], 0, 1, dt1)

	#initial extinction probability at 
	# for second integration under MEA type "recompute"

	e1 <- E_func(pars0["lambda"], pars0["mu"], 0, dt1)

	# data probability for secont (x0) segment:
	d2 <- D_func(pars0["lambda"], pars0["mu"], e1, d1, dt0)

	# extinction at the root
	e2 <- E_func(pars0["lambda"], pars0["mu"], 0, total_age)

	# probability:
	return(as.vector(d2 / (1 - e2)))
	
} 














