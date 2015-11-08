# This script illustrates issues with combining extinction probabilities at nodes.
# 
#
#
#

# likelihood functions for the constant-rate birth-death process
#    used to compute speciation and extinction probabilities
#    on individual branches, given initial E and D probabilities.

# Probability of lineage extinction before the present:
E_func <- function(lam, mu, E0, dt) {
	
	num <- (1 - E0) * (lam - mu);
	denom <- (1 - E0) * lam - (mu - lam * E0) * exp(-(lam - mu) * dt);
	
	return( 1 - num / denom)	
 
}

# Probability of the data.
D_func <- function(lam, mu, E0, D0, dt) {
	
	r <- lam - mu
	num <- (D0 * r^2) * exp((-r) * dt);
	denom <- ( (lam - lam * E0 + exp((-r) * dt) * (lam * E0 - mu)  )  ) ^ 2;
	return(num / denom) 
 
}


# EXERCISE 1: LIKELIHOOD OF SIMPLE 2 TAXON TREE WITH STEM BRANCH
# AND RATE SHIFT
# Compute likelihood of a simple tree with one internal node:
# (A:99,B:99):1
# 
#  Imagine shift happens IMMEDIATELY
#  At origin of lineage A. 
#

# What is likelihood of tree as a function of mu_root?
# Hold lambda_root = 0.5
#
#
mu_vec <- seq(0.0001, 1, length.out=100)
lmat <- matrix(NA, nrow=100, ncol=2)

# A branch likelihood is 0


for (i in 1:length(mu_vec)){
	
	
	lhA <- log(D_func(0.5, mu_vec[i], 0, 1, 99))
	 
	exprob_A <- 0
	exprob_B <- E_func(0.5, mu_vec[i], 0, 99)
	
	exprob_AB_1 <- exprob_A * exprob_B
	exprob_AB_2 <- exprob_B
	
	exprob_root_1 <- E_func(0.5, mu_vec[i], exprob_AB_1, 1)
	exprob_root_2 <- E_func(0.5, mu_vec[i], exprob_AB_2, 1)
		
 	D1 <- log(D_func(0.5, mu_vec[i], 0, 1, 99)) 
 	D1 <- D1 + log(D_func(0.5, mu_vec[i], exprob_AB_1, 1, 1))
 	D1 <- D1 + log(0.5) # likelihood of the speciation event
 	
 	D2 <- log(D_func(0.5, mu_vec[i], 0, 1, 100))
 	D2 <- D2 + log(0.5)	
	
	lmat[i,1] <- D1 - log(1 - exprob_root_1)
	lmat[i,2] <- D2 - log(1 - exprob_root_2)
	
}


png(height = 600, width=600, file = "likelihood_nodecombine.png")
plot.new()
par(mar=c(6,6,1,1))
plot.window(xlim=c(0, 1), ylim=c(-50, 15))
lines(mu_vec, lmat[,1], lwd=5, col="red")
lines(mu_vec, lmat[,2], lwd= 5, col="blue")
arrows(x0 = 0.6, x1= 0.85, y0 = 10, y1 = -0.5, lwd=2, length=0.15, col="blue")

arrows(x0 = 0.5, x1 = 0.5, y0=-40, y1=-51, lwd=2, length=0.15)

axis(1, at=seq(-0.2, 1.2, by =0.2), cex.axis=1.3)
axis(2, at=seq(-60,10, by=10), las=1, cex.axis=1.3)
text(x=0.5, y=-38, label="Root speciation rate", cex=1.7, font=3)

mtext("Extinction rate for root process", 1, cex=1.9, line=4)
mtext("log-likelihood of tree", 2, cex=1.9, line=4)
text(x=0.6, y=10, label="Numerical failure", pos=2, cex=1.7, font=3, col="blue")
 


dev.off()


# EXERCISE 2: LIKELIHOOD OF 4-TAXON TREE AS SHOWN ON WEBSITE
#
# this pertains to a tree of the form
# ((A:99,B:99):1,(C:99,D:99):1);
# We will assume that shifts happen on branches A and C
# and that they happen at exactly 1 time unit into the branch.
# thus, on the terminal branch A and C, the lineage is diversifying under
# root process parameters for 1.0 time unit, and the new parameters for 98.0 time units.
# The internal branches are under the parent process
# This is the 4-taxon tree that is illustrated on the website linking this code.

# We will assume that the parameters on each of the shift branches are identical. Thus, you have 
# for the constant-rate birth-death process, 4 parameters:
# lambda_root, mu_root, lambda_shift, mu_shift.

# Here we will create two likelihood functions,
# hard-coding the tree structure, but accepting a 4-parameter vector 
# of the 4 parameters in the order:
# lambda_root, mu_root, lambda_shift, mu_shift.

# The likelihood functions will differ in how they combine probabilities at nodes.

# This likelihood function returns the likelihood and can be toggled 
# between 2 methods of handling the extinction probabilities 
# at nodes
# combine.type = favor_parent 
#							will always favor the parent process
# combine.type = multiply 
#							will multiply the values E(t) from right 
#							and left branches
#
#

likelihood4taxon <- function(p, combine.type){
 
	
	lam_r <- p[1]
	mu_r <- p[2]
	lam_s <- p[3]
	mu_s <- p[4]
	
	# note that branches B and D are mirror images, 
	# as are A and C. Thus we will simply compute for A and B
	# and double them.
	
	# probs for B branch:
	d_b <- D_func(lam_r, mu_r, 0, 1, 99)
	e_b <- E_func(lam_r, mu_r, 0, 99)
	
	# probs for A branch:
	d_t1 <- D_func(lam_s, mu_s, 0, 1, 98)
	e_t1 <- E_func(lam_s, mu_s, 0, 98)
	d_a <- D_func(lam_r, mu_r, e_t1, d_t1, 1)
	e_a <- E_func(lam_r, mu_r, e_t1, 1)
	
	
 	e_basal <- 0
 	if (combine.type == "favor_parent"){
 		e_basal <- e_b
 	}else if (combine.type == "multiply"){
 		e_basal <- e_a * e_b
 	}else{
 		stop("invalid combine.type")
 	}
	
	d_basal <- D_func(lam_r, mu_r, e_basal, 1, 1)
	e_basal <- E_func(lam_r, mu_r, e_basal, 1)
	
	
	loglik <- 2 * log(d_b) + 2 * log(d_a) + 2 * log(d_basal)
	loglik <- loglik + 2*log(lam_r)
	
	# adding constant term for comparison with diversitree:
	loglik <- loglik + sum(log(2:3))
	
	# conditioning on survival:
	loglik <- loglik - 2 * log(1 - e_basal)
	return(loglik)
}


### LIKELIHOOD UNDER THE 2-PARAMETER MODEL: NO SHIFTS

### Here we will fit a 2 parameter model to the tree, just assuming a 
# constant rate birth-death process.
# We will then compare to a 4-parameter process that has the shift parameters
# on lineage A and lineage C.

lwrap2 <- function(pars, combine.type){
	# optimizing on logscale...
	
	pars <- exp(c(pars, pars))
	return(likelihood4taxon(pars, combine.type))
}

# This is the ML estimate for parameters
# under the 2 parameter model. Note that this does not have a shift,
# so we must use the favor_parent option for E(t) handling at nodes.
res2 <- optim(log(c(1,0.5)), lwrap2, method="Nelder", control=list(fnscale = -1), combine.type = "favor_parent")

# likelihood at the maximum: 
res2$value

# the ML params
pars2 <- exp(res2$par)
names(pars2) <- c("lambda", "mu")

# Just to check, compare to diversitree
# which can also compute the likelihood of the constant-rate birth-death process

tree <- read.tree(text = "((A:99,B:99):1,(C:99,D:99):1);")
library(diversitree)
lfx <- make.bd(tree)

# diversitree likelihood includes a constant that 
lfx(pars2)

lfx(pars2) - res2$value

# a very small number. So Diversitree 
# computes identical likelihood for the constant-rate process. 

### LIKELIHOOD UNDER THE 4-PARAMETER MODEL: RATE SHIFTS + E(T) TAKEN FOR PARENT PROCESS
#
# Define a new wrapper

lwrap4 <- function(pars, combine.type){
	# optimizing on logscale...
	
	pars <- exp(pars)
	return(likelihood4taxon(pars, combine.type))
}

init <- c(0.1,0.01, 0.1,0.01)

# This optimization is going to fail often due to numerical issues
# when the extinction probability becomes too close to 1
# so we will do many optimizations and save results.
# We will do w random starting parameters,
# including many where the initial extinction rate is higher than 
# speciation

N_OPTIMIZE <- 100

xx_fp <- matrix(NA, nrow= N_OPTIMIZE, ncol=5) # fp for "favor_parent"
colnames(xx_fp) <- c("logLik", "lambda_r", "mu_r", "lambda_s", "mu_s")

for (i in 1:N_OPTIMIZE){
	
	init <- runif(4)		

	tmp <- try(optim(log(init), lwrap4, method="Nelder", combine.type = "favor_parent", control=list(fnscale = -1)))
	if (class(tmp) != "try-error"){
		xx_fp[i,1] <- tmp$value
		xx_fp[i,2:5] <- exp(tmp$par)
	}
	
}

# get rid of NAs...
xx_fp <- xx_fp[!is.na(xx_fp[,1]), ]

best <- which(xx_fp[,1] == max(xx_fp[,1]))[1]

ml_fp <- xx_fp[best,]

### Look at this parameter set and what appears to be wrong:
### HUGE improvement in likelihood for this model:
### Compare 2 parameter vs 4 parameter model:
### For a tiny 4 taxon tree 

ml_fp[1] - res2$value

## But it is especially illuminating to consider the parameters:
## With the "favor_parent" option, we estimate that most of the tree
## diversified under a process with much higher extinction than speciation
## specifically, the basal branches + lineages B and D diversified with relative
# extinction rate of:

ml_fp[3] / ml_fp[2]

# So the ML estimate here is that extinction was 150% the speciation rate. 
# Very strange.
#
# For comparison, you can ask: just how likely would such a tree be 
# under a single constant-rate birth-death process with those parameters:

lfx(ml_fp[2:3])

# Not very likely. 
# Diversitree is giving a loglikelihood difference of about 70 for this.
# Let us compare to what we get if we combine extinction at nodes 
# by multiplying, as is now the default in BAMM:
#
#  
# THIS VALUE (FROM FAVOR_PARENT) MAY BE THEORETICALLY IMPOSSIBLE,
#
#
# If the Diversitree log-likelihood is -70 for these parameters in this 
# 4 taxon tree under a 2 parameter birth-death process
#
# then the theoretical maximum improvement to the log-likelihood
# under these same parameters, with rate shifts on branches A and C 
# is zero 
#
# in the limit as t -> infinity and lambda -> 0, 
# the probability of a single branch reaches a maximum of 1.
# Thus, even after you add rate shifts to the branches for A and C,
# the best you can do, for this set of parameters, should only be -70.
# But somehow the model is giving us a value of -2. 
# 

### LIKELIHOOD UNDER THE 4-PARAMETER MODEL: RATE SHIFTS + E(T) multiplied together
#

xx_co <- matrix(NA, nrow= N_OPTIMIZE, ncol=5) # co for "combine"
colnames(xx_co) <- c("logLik", "lambda_r", "mu_r", "lambda_s", "mu_s")

for (i in 1:N_OPTIMIZE){
	
	init <- runif(4)		

	tmp <- try(optim(log(init), lwrap4, method="Nelder", combine.type = "multiply", control=list(fnscale = -1)))
	if (class(tmp) != "try-error"){
		xx_co[i,1] <- tmp$value
		xx_co[i,2:5] <- exp(tmp$par)
	}
	
}


xx_co <- xx_co[!is.na(xx_co[,1]), ]

best <- which(xx_co[,1] == max(xx_co[,1]))[1]

ml_co <- xx_co[best,]

# And now the log-likelihood is just -9.43
# The Diversitree log-likelihood, 
# if the root process applied to the entire tree
# is this:

lfx(ml_co[2:3])

# I get -11.38.
# So, allowing separate shifts on A and C gives a modest improvement 
# here, to -9.43.

 
AIC_2 <- as.numeric(-2 * res2$value + 4)
AIC_4_fp <- as.numeric(-2 * ml_fp[1] + 8)
AIC_4_co <- as.numeric(-2 * ml_co[1] + 8)

# No improvement in AIC for slapping on shifts to terminal branches 
# if multiplying extinction probs
AIC_2  - AIC_4_co 

# Big improvement if you use the favor_parent approach
AIC_2  - AIC_4_fp 



#################### 
 
 
 
 




















