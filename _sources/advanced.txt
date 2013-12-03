.. highlight:: none

Advanced analysis options
============================

Modeling *less complex* evolutionary scenarios
------------------------------

The basic BAMM model is fairly complex, as it allows rate variation through time and among lineages. However, you can easily modify the control file for BAMM to perform several simpler types of analyses.

Constant rate birth-death model
************************
To perform Bayesian inference on your data under a constant-rate birth-death process, you can simply tell BAMM not to perform any MCMC update moves that are not part of the constant-rate birth-death process. Specifically, you should (i) not add rate regimes to the tree, and (ii) not update the parameter controlling speciation rate variation through time. You can do this by making sure the MCMC move frequencies are set as follows::

	updateRateEventNumber = 0
	updateRateEventPosition = 0
	updateRateEventRate = 0
	updateRateLambda0 = 1
	updateRateLambdaShift = 0
	updateRateMu0 = 1
 

Pure-birth model
************************
To run a pure-birth only model, with no extinction, you just turn off the extinction update::
	
	updateRateMu0 = 0
	
However, you must be careful to ensure that the initial value of extinction is set to zero. Since you are no longer updating this parameter through MCMC, whatever value it starts with is the value it will be stuck with::
	
	muInit0 = 0.0
	 
Note that you cannot account for incomplete sampling under a pure-birth model (modeling incomplete sampling is mathematically identical to allowing a particular type of extinction).

MEDUSA-like model
************************
To run a Bayesian MEDUSA-like model, where the rate of speciation and extinction is constant within specific shift regimes, you can set your MCMC move frequencies as follows::
	
	updateRateEventNumber = 0.1
	updateRateEventPosition = 1
	updateRateEventRate = 1
	updateRateLambda0 = 1
	updateRateLambdaShift = 0
	updateRateMu0 = 1
	
And of course, we need to ensure that::

	lambdaShift0 = 0	
	
Here, we are simply setting the time-variation parameter (lambdaShift) of the BAMM model to zero, and also setting the update frequency for that parameter to zero. 

BAMM is sufficiently flexible as to allow a number of permutations on these general themes. In addition, the modifications to model setup described here also apply to trait evolution, where you could just as easily constrain a phenotypic analysis to involve only time-invariant Brownian motion processes (similar to the *Auteur* package for R), with the following code::
	
	updateRateBetaShift = 0.0
 	
 	betaShiftInit = 0.0
 
 
 
Accounting for phylogenetic uncertainty
------------------------------

Some researchers consider it important to account for phylogenetic uncertainty when performing macroevolutionary analyses. At present, there is no direct way of accounting for phylogenetic uncertainty in BAMM itself. It remains unclear whether phylogenetic uncertainty generally matters for the sorts of conclusions obtained with BAMM. My (DLR) personal view is that phylogenetic uncertainty is very much an issue for **some types** of results obtained using BAMM (and other programs), and (usually) not an issue at all for many other types of results. 

When does phylogenetic uncertainty **not** matter? For general inference on the overall tempo and mode of diversification, it is quite unlikely that - in general - your focal tree (say, MCC tree from BEAST, or ML tree from RAxML) is unlikely to be *so bad* that your broad-scale inferences in evolutionary dynamics are inaccurate. One reason for this is the somewhat paradoxical observation that confidence in *macroevolutionary conclusions* can be negatively correlated with confidence in *phylogenetic conclusions*. For example, consider a phylogeny showing a pattern of an early burst in lineage diversification, such that most major lineages arose during a brief period of time. You might never be able to resolve the *precise* order of branching of those lineages, and as such, you will always have a tree that is poorly resolved at the base. However, you might nonetheless be extremely confident that branch lengths are short near the base of the radiation (indeed, this is why you have low confidence in your topology!), and this means that your inferences on speciation rates themselves might be extremely robust. If you consider speciation in whales, as shown :ref:`here<whalefig>` and :ref:`here<rttwhale>`, phylogenetic uncertainty isn't going to change the big-picture conclusions: there was clearly a massive spike in evolutionary rates in some ancestral lineage leading to, or immediately nested within, the dolphin clade.

Phylogenetic uncertainty will matter if you do in fact care about *specific* aspects of changes in evolutionary dynamics. If you really care about the *precise* location of the shift in evolutionary dynamics, then the exact sequence of branching at the base of the dolphin radiation (to continue with the aforementioned example) **will** matter. Please keep in mind, however, that the BAMM model (and all other models), are merely statistical models that have imposed on the data. So, excessively fretting about whether the true shift in evolutionary dynamics occurred on branch *A* or branch *B* is somewhat unproductive, because the notion of a discrete shift is itself an assumption of the model we are using for inference.

Although BAMM does not directly allow modeling phylogenetic uncertainty, it is straightforward to perform BAMM analyses across distributions of phylogenies taken from a Bayesian analysis. Below, we demonstrate how you can use your bash shell (on the OSX or Linux operating systems) to perform a BAMM analysis across a sample of trees.::

	Add bash script stuff here, maybe a link to ppss as well.
	
	
	
	




 