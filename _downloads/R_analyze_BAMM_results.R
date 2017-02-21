
# This script assumes that you have run BAMM
#   on trees simulated with MEA's tree simulator
#   with random starting seeds of 8, 9, and 15
#
#  Why these particular seeds?
# 
#  In the first 20 integer-valued seeds (1 to 20)
#   these 3 are the only values that gave tree sizes 
#   that were close to the size selection threshold imposed 
#   by MEA
#   Most other seeds generated trees that had 
#   more than 1000 taxa 
# 
#  This script assumes that you have:
#   1) Generated a tree with MEAs simulator and a starting seed
#   2) analyzed this tree with BAMM
#   3) Understand how to specify the filenames to access the output from BAMM
#   4) Have generated an "event data" file for the TRUE rate shift data
#          from the MEA simulator. 
#          In the script file "R_simulate_tree_from_seed.R", the 
#          true event data is written to file using a naming convention
#          sX_true_eventdata.txt, where X is the target random seed.
#
# ----------------------------------------------------------------
 
library(ape)
library(BAMMtools)
library(coda)
source("R_test_BAMM_functions.R")

# ----------------------------------------------------------------
# SEED 8
# A relatively "high-power" tree  

# Read in tree, mcmc data, and BAMM output
tree8    <- read.tree("tree_s8.tre")
mcmc8    <- read.csv( "s8_mcmc_out.txt" , stringsAsFactors=F)
ed8      <- getEventData(tree8,  "s8_event_data.txt", burnin=0.1, nsamples=200)
true_ed8 <- getEventData(tree8, "s8_true_eventdata.txt")

# Now we count the number of tips in each shift regime

ss8 <- shiftRegimeData(true_ed8)

# look at ss:
ss8

# This is all the shift data from the simulated tree, including 
#  true parameters, time and location of shift, and the 
#  number of tips in each shift regime.


# There are 69 tips in the root regime
# 52 tips in the 2nd largest shift regime
# 9 tips in the 3rd largest regime (row 5)
# This suggests we might have *some* power to infer rate variation
 
# ----------------------------------------------------------------
# Make a pair of plots for visually inspecting BAMM's performance
 
# quartz.options(height=7, width=11) # optional OSX plot setup command 
plot.new()
par(mfrow=c(1, 2))
BAMMplot <- plot.bammdata(true_ed8, spex="s", breaksmethod="linear", lwd=2, tau=0.003)
mtext("True rates", side=3, cex=1.5)
addBAMMlegend(BAMMplot, location="left")
plot.bammdata(ed8, colorbreaks=BAMMplot$colorbreaks, spex="s", lwd=2, tau=0.003)
mtext("BAMM estimated rates", side=3, cex=1.5)
 
 
# ----------------------------------------------------------------
# Explicit analysis of how well BAMM reconstructed 
#   "branch specific" speciation rates:
 
# First, we compute a tree where the branch lengths
#   are equal to the BAMM-estimated speciation rates,
#   using a function from BAMMtools

rates8   <- getMeanBranchLengthTree(ed8)$phy

# Here is the corresponding tree where branch lengths equal TRUE rates
rates8true <- getMeanBranchLengthTree(true_ed8)$phy

# First, we fit a linear model to the estimated vs true 
#   branch-specific rates
fit <- lm(rates8$edge.length ~ rates8true$edge.length)


# Now, we plot the rates:
# First, we set up a simple plotting function to make all of this easier
# This fx takes a single argument: the maximum rate value
# eg max of x or y axis:

plotSetup <- function(mx = 0.8){
	plot.new()
	par(oma=c(1,1,1,1))
	par(mar=c(5,5,1,1))
	plot.window(xlim=c(0, mx), ylim=c(0, mx), asp=1)	
	abline(0, 1, lwd=3, col="gray50")
	axis(1, at=seq(-0.2, 1.2, by=0.2))
	axis(2, at=seq(-0.2, 1.2, by=0.2), las=1)
	mtext("True rate", side=1, line=3, cex=1.4)
	mtext("BAMM estimated rate", side = 2, line = 3, cex = 1.4)
	
	lines(x=c(0,0.1), y=c(0.6, 0.6), lwd=3, col="gray50")
	text(x=0.1, y=0.6, pos=4, label = "1:1 line", font = 3)
	lines(x=c(0,0.1), y=c(0.55, 0.55), lty="dotted", lwd=3)
	text(x=0.1, y=0.55, pos=4, label = "BAMM fit", font=3)
		
	
}

# There is a huge amount of 
#  overplotting in the following plots, 
#   since most taxa have identical rates.
#  So we will "jitter" the points slightly so you can see
#  roughly how many branches fall out together

j <- function(x){
	return(jitter(x, amount = 0.015))
}
 

# quartz.options(height = 6, width=6)
plotSetup()

# add fitted line:
abline(fit$coefficients[1], fit$coefficients[2], lwd=3, lty="dotted")
 
# plotting points, using jitter function from above: 
#points(j(rates8true$edge.length), j(rates8$edge.length), pch=21, bg="coral") 


colset <- colorByRegime(true_ed8)
# plotting points, using jitter function from above: 
points(j(rates8true$edge.length), j(rates8$edge.length), pch=21, bg=colset) 

# ----------------------------------------------------------------
# Average branch rates for nodes 148 and 240:
mean(extract.clade(as.phylo(rates8), node = 148)$edge.length)
# gives 0.56

mean(extract.clade(as.phylo(rates8), node = 240)$edge.length)
# gives 0.52

# ----------------------------------------------------------------
# Diversitree rates for clades from nodes 148 and 240:
birthdeath_fit(tree8, node = 148) 

# Diversitree rates for clades from nodes 148 and 240:
birthdeath_fit(tree8, node = 240) 
  
 
  
 
 # ----------------------------------------------------------------
# Rates by clade for clades from nodes 148 and 240:

cr148 <- getCladeRates(ed8, node = 148)
mean(cr148$lambda)
# BAMM estimate = 0.46 vs 0.39 true

cr240 <- getCladeRates(ed8, node = 240)
mean(cr240$lambda)
# BAMM estimate = 0.43 vs 0.41 true

#
# ----------------------------------------------------------------
# SEED 9
# A low-power tree (would have been included in MEA)

# Read in tree, mcmc data, and BAMM output
tree9    <- read.tree("tree_s9.tre")
mcmc9    <- read.csv( "s9_mcmc_out.txt" , stringsAsFactors=F)
ed9      <- getEventData(tree9,  "s9_event_data.txt", burnin=0.1, nsamples=200)
true_ed9 <- getEventData(tree9, "s9_true_eventdata.txt")

# Now we count the number of tips in each shift regime

ss9 <- shiftRegimeData(true_ed9)

# look at ss:
ss9

# This is all the shift data from the simulated tree, including 
#  true parameters, time and location of shift, and the 
#  number of tips in each shift regime.


# There are 81 tips in the root regime
# 3 tips in the 2nd largest shift regime
# 1 tip each of 2 other shift regimes.
# It would be surprising if BAMM could detect shifts on this tree:
 
# ----------------------------------------------------------------
# Make a pair of plots for visually inspecting BAMM's performance
 
# quartz.options(height=7, width=11) # optional OSX plot setup command 
plot.new()
par(mfrow=c(1, 2))
BAMMplot <- plot.bammdata(true_ed9, spex="s", breaksmethod="linear", lwd=2.4, tau=0.003)
addBAMMshifts(true_ed9, par.reset=F, cex=2, bg = "black")
mtext("True rates", side=3, cex=1.5)
addBAMMlegend(BAMMplot, location="left")
plot.bammdata(ed9, colorbreaks=BAMMplot$colorbreaks, spex="s", lwd=2, tau=0.003)
mtext("BAMM estimated rates", side=3, cex=1.5)

# ----------------------------------------------------------------
# mean whole-tree rates

mean(getCladeRates(ed9)$lambda)
# I get 0.20

# the "true" mean:
mean(getCladeRates(true_ed9)$lambda)
# I get 0.21


# ----------------------------------------------------------------
# Explicit analysis of how well BAMM reconstructed 
#   "branch specific" speciation rates:
 
# First, we compute a tree where the branch lengths
#   are equal to the BAMM-estimated speciation rates,
#   using a function from BAMMtools

rates9   <- getMeanBranchLengthTree(ed9)$phy

# Here is the corresponding tree where branch lengths equal TRUE rates
rates9true <- getMeanBranchLengthTree(true_ed9)$phy

# Fit a linear model to the estimated vs true 
#   branch-specific rates
fit9 <- lm(rates9$edge.length ~ rates9true$edge.length)

summary(fit9)

# quartz.options(height = 6, width=6)
plotSetup()

# add fitted line:
abline(fit9$coefficients[1], fit9$coefficients[2], lwd=3, lty="dotted")
 
colset <- colorByRegime(true_ed9)
# plotting points, using jitter function from above: 
points(j(rates9true$edge.length), j(rates9$edge.length), pch=21, bg=colset) 
 
 


# ----------------------------------------------------------------
# ----------------------------------------------------------------
# SEED 15
# A high-power tree (excluded from MEA)

# Read in tree, mcmc data, and BAMM output
tree15    <- read.tree("tree_s15.tre")
mcmc15   <- read.csv( "s15_mcmc_out.txt" , stringsAsFactors=F)
ed15      <- getEventData(tree15,  "s15_event_data.txt", burnin=0.1, nsamples=200)
true_ed15 <- getEventData(tree15, "s15_true_eventdata.txt")

# Now we count the number of tips in each shift regime

ss15 <- shiftRegimeData(true_ed15)


# ----------------------------------------------------------------
# Make a pair of plots for visually inspecting BAMM's performance
 
# quartz.options(height=7, width=11) # optional OSX plot setup command 
plot.new()
par(mfrow=c(1, 2))
BAMMplot <- plot.bammdata(true_ed15, spex="s", breaksmethod="linear", lwd=2, tau=0.003)
mtext("True rates", side=3, cex=1.5)
addBAMMlegend(BAMMplot, location="left")
plot.bammdata(ed15, colorbreaks=BAMMplot$colorbreaks, spex="s", lwd=2, tau=0.003)
mtext("BAMM estimated rates", side=3, cex=1.5)


rates15   <- getMeanBranchLengthTree(ed15)$phy

# Here is the corresponding tree where branch lengths equal TRUE rates
rates15true <- getMeanBranchLengthTree(true_ed15)$phy

# Fit a linear model to the estimated vs true 
#   branch-specific rates
fit <- lm(rates15$edge.length ~ rates15true$edge.length)

# quartz.options(height = 6, width=6)
plotSetup(1)

# add fitted line:
abline(fit$coefficients[1], fit$coefficients[2], lwd=3, lty="dotted")
 
colset <- colorByRegime(true_ed15)
# plotting points, using jitter function from above: 
points(j(rates15true$edge.length), j(rates15$edge.length), pch=21, bg=colset) 
 
# ----------------------------------------------------------------
# Is BAMM doing good or bad? Let's compare to diversitree 
#   on these same sets of shift clades:
 
# Node that there are essentially 3 shift regimes on this tree:
# the root regime:
ss15[1, ]

# This regime with 102 tips and fast speciation:
ss15[3, ]

# and this regime with even faster speciaion and 39 tips:
ss15[4, ]

# You can see from the plot above that BAMM overestimated one shift clade
#   and underestimated the other. 

# regime 2: true lambda = 0.49, mu = 0.08
# BAMM estimate for this clade (note that this is node 411)

cr <- getCladeRates(ed15, node = 411)
mean(cr$lambda)   
mean(cr$mu)

# BAMM estimates are almost spot-on


# how about the other regime:
ss[4,]
# here, lambda = 0.60, mu = 0; node = 345
cr <- getCladeRates(ed15, node = 345)
mean(cr$lambda)   # I get 0.8
mean(cr$mu)       # I get 0.22

# BAMM overestimated speciation and extinction, but here the 
#   net diversification rate is very accurate (0.58 vs 0.6)

# How well does Diversitree do for these 2 shift clades?
# You can use the birthdeath_fit function that we have created 
#  for this. It is just a wrapper for diversitree likelihood functions 
#  for the constant-rate birth-death process, but lets you extract and fit
#   the model to a specific subclade

birthdeath_fit(tree15, node = 411)
# I get 0.451 and 0, which is close to the BAMM estimate
# and slightly further away (but still close) to true value

birthdeath_fit(tree15, node = 345)
# For the second shift clade:
# I get 0.81 and 0.08, The speciation estimate is almost identical to BAMM
#   and overestimates the true value, 
#   but extinction is closer.
#   However, BAMM's net diversification rate estimate is more accurate in this case.
 







