

#############################################################
#
#	getEventData(....)
#
#	eventfilename		=	file where event data are stored, in bamm output format
#	nsamples			=	number of samples from posterior to include. 
#							if NULL, includes ALL samples (could take long time to run)
#	phy					=	The model tree (time calibrated tree analyzed w bamm)
#	verbose				=	Verbose print output for bug tracking
#	burnin				=	How many samples from posterior to discard
#							Uses fraction (e.g, 0.25 = 25% discarded)
#	type				=	specifies whether eventfilename refers to trait or diversification data
#	header				=	Boolean to flag whether eventfilename contains a header
getEventData <- function(phy, eventfilename, burnin=0, nsamples = NULL, verbose=FALSE, assign.type = 'new_way', type = 'diversification', header=TRUE){
	
	if (type != 'diversification' & type != 'traits'){
		stop("Invalid 'type' specification. Should be 'diversification' or 'traits'");
	}
	
	if (any(is.null(c(phy$begin, phy$end)))){
		phy <- getStartStopTimes(phy);
	}
		
	
	bt <- branching.times(phy);
	
	eventVectors <- list();
	eventData <- list();
	tipStates <- list();
	eventBranchSegs <- list();
	
	tipLambda 	<- list();
 
	cat("Reading event datafile: ", eventfilename, "\n\t\t...........");
 	x <- read.csv(eventfilename, header=header, stringsAsFactors=FALSE);
 	uniquegens <- sort(unique(x[,1]));
 	cat("\nRead a total of ", length(uniquegens), " samples from posterior\n");
 	
 	samplestart <- uniquegens[floor(burnin*length(uniquegens))];
 	if(!length(samplestart))
 	{
 		samplestart <- 0;
 	}
 	uniquegens <- uniquegens[uniquegens >= samplestart];
 
 	if (is.null(nsamples)){
 		nsamples <- length(uniquegens);
 	}else if (nsamples > length(uniquegens)){
 		nsamples <- length(uniquegens);
 	}
	
	goodsamples <- uniquegens[seq.int(1, length(uniquegens), length.out=nsamples)];
 	 
 	cat('\nDiscarded as burnin: GENERATIONS < ', goodsamples[1]);
 	cat("\nAnalyzing ", length(goodsamples), " samples from posterior\n");

 	numberEvents <- length(goodsamples); # vector to hold number of events per sample
 
 	cat('\nSetting recursive sequence on tree...\n');
 	
 	if (assign.type=='new_way')
 		phy <- getRecursiveSequence(phy);
 	cat('\nDone with recursive sequence\n\n');
 
	meanTipMu <- numeric(length(phy$tip.label));
	
 	meanTipLambda <- numeric(length(phy$tip.label)); 
 
 	for (i in 1:length(goodsamples)){
  		
  		tmpEvents <- x[x[,1] == uniquegens[i], ];
		
		if (verbose)
			cat('Processing event: ', i, '\n');		
 
		t1 <- tmpEvents[,2]; # desc 1
 		t2 <- tmpEvents[,3]; # desc 2
 		tm <- tmpEvents[,4]; # abs time of event
 		lam1 <- tmpEvents[,5]; # lambda parameter 1
 		lam2 <- tmpEvents[,6]; # lambda parameter 2
 		if(type == 'diversification'){
			mu1 <- tmpEvents[, 7]; # mu parameter 1
 			mu2 <- tmpEvents[, 8]; #mu parameter 2 
		}
		else{ #for bamm trait data we set the mu columns to zero because those params don't exist
			mu1 <- rep(0, nrow(tmpEvents)); 
 			mu2 <- rep(0, nrow(tmpEvents)); 
		}
		tipMu <- list();	
		
 		
 		# Get subtending node for each event:
 		nodeVec <- numeric(nrow(tmpEvents));
 		for (k in 1:length(nodeVec)){
 			if (is.na(t2[k])){
 				# Node is a tip
 				nodeVec[k] <- which(phy$tip.label == t1[k]);	
 			}else{
 				tipnodes <- c(which(phy$tip.label == t1[k]), which(phy$tip.label == t2[k]));
 				nodeVec[k] <- getMRCA(phy, tipnodes);
 			}
 		}
 		
 		if (sum(nodeVec == 0)  > 0){
			stop('Failed to assign event to node\n');
		}
		
		# make a dataframe:
		dftemp <- data.frame(node=nodeVec, time=tm, lam1=lam1, lam2=lam2, mu1=mu1, mu2=mu2, stringsAsFactors=FALSE);
		
		
		dftemp <- dftemp[order(dftemp$time), ];
		dftemp$index <- 1:nrow(dftemp);
		rownames(dftemp) <- NULL;
		
		tphy <- phy;
		#tphy$statevec <- numeric(nrow(phy$edge));	
		tphy$statevec <- rep(1, nrow(tphy$edge));

		if (assign.type == 'old_way'){
			
			for (k in 1:nrow(dftemp)){
				tphy <- recursivelySetNodeStates(tphy, node=dftemp$node[k], state=k);
			}			
			
		}else if (assign.type == 'new_way'){
			if (nrow(dftemp) > 1){
				for (k in 2:nrow(dftemp)){
					
					s1 <- which(tphy$downseq == dftemp$node[k]);
					s2 <- which(tphy$downseq == tphy$lastvisit[dftemp$node[k]]);
					descSet <- tphy$downseq[s1:s2];
					isDescendantNode <- tphy$edge[,2] %in% descSet;				
					tphy$statevec[isDescendantNode] <- k;
				}
							
			}

		}

 		tmpEventSegMat <- matrix(0, nrow=(max(phy$edge) + nrow(dftemp) - 2), ncol=4);
		non.root <- c(1:length(phy$tip.label), (length(phy$tip.label)+2):max(phy$edge));
		pos <- 1;	
		
		for (k in non.root){
 			events.on.branch <- dftemp[dftemp$node == k, ];
			if (nrow(events.on.branch) == 0){
				tmpEventSegMat[pos,1] <- k;
				tmpEventSegMat[pos,2] <- tphy$begin[tphy$edge[,2] == k];
				tmpEventSegMat[pos,3] <- tphy$end[tphy$edge[,2] == k];
				tmpEventSegMat[pos,4] <- tphy$statevec[tphy$edge[,2] == k];
				
			}else{
				events.on.branch <- events.on.branch[order(events.on.branch$time), ];
				
				fBranch <- phy$edge[,2] == k;
 				start.times <- c(phy$begin[fBranch], events.on.branch$time);
				stop.times <- c(events.on.branch$time, phy$end[fBranch]);
				parent <- phy$edge[,1][phy$edge[,2] == k];
				if (parent == (length(phy$tip.label) + 1)){
					# Parent is root:
					proc.set <- c(1, events.on.branch$index);	
				}else{
					proc.set <- c(tphy$statevec[tphy$edge[,2] == parent], events.on.branch$index);			
				}
				
 				zzindex <- pos:(pos+nrow(events.on.branch));	
				
				tmpEventSegMat[zzindex, 1] <- rep(k, length(zzindex));
				tmpEventSegMat[zzindex, 2] <- start.times;
				tmpEventSegMat[zzindex, 3] <- stop.times;
				tmpEventSegMat[zzindex, 4] <- proc.set;
				
			}
			
			pos <- pos + 1 + nrow(events.on.branch);
		} 	
 	
 		eventBranchSegs[[i]] <- tmpEventSegMat;

		tipstates <- numeric(length(phy$tip.label));
		for (k in 1:length(tipstates)){
			tipstates[k] <- tphy$statevec[tphy$edge[,2] == k];
		}
 		
 		### Compute tip rates:
 
		stoptime <- max(branching.times(phy));
		
		tiplam <- dftemp$lam1[tipstates] * exp(dftemp$lam2[tipstates] * (stoptime - dftemp$time[tipstates]));
		tipmu <- dftemp$mu1[tipstates];
		
		meanTipMu <- meanTipMu + tipmu/nsamples;
		meanTipLambda <- meanTipLambda + tiplam/nsamples;
		
		
		
		### List assignments and metadata across all events:
		eventData[[i]] <- dftemp;	
		eventVectors[[i]]  <- tphy$statevec;
		numberEvents[i] <- nrow(dftemp);
		tipStates[[i]] <- tipstates;
		
		tipLambda[[i]] <- tiplam;
		tipMu[[i]] <- tipmu;	
 	}
 	
	phy$numberEvents <- numberEvents;
	phy$eventData <- eventData;
	phy$eventVectors <- eventVectors;
	phy$tipStates <- tipStates;
	phy$tipLambda <- tipLambda;
	phy$meanTipLambda <- meanTipLambda;
	phy$eventBranchSegs <- eventBranchSegs; 	
	phy$tipMu <- tipMu;
	phy$meanTipMu <- meanTipMu;
	if(type == 'diversification'){	
		phy$type = 'diversification';
	}
	else{
		phy$type = 'traits';	
	}
 	
	# Inherits attributes of class phylo
	# plus adds new class: bamm-data
	class(phy) <- c('phylo', 'bamm-data');
	return(phy);
}

#############################################################
#
#	branchMeanRateExponential(....)
#
#	Computes time-averaged rate of eponential process
#	Vectorized


branchMeanRateExponential <- function(t1, t2, p1, p2){
	if (p2 == 0){
		return(p1);
	}else{
		(p1/p2)*(exp(p2*t2) - exp(p2*t1)) / (t2 - t1);
	}
}
branchMeanRateExponential <- Vectorize(branchMeanRateExponential);

#############################################################
#
#	exponentialRate(....)
#
#	Gets point estimate of evolutionary rate
#	Vectorized

exponentialRate <- function(t1, p1, p2){
	(p1 * exp(p2 * t1));
}
exponentialRate <- Vectorize(exponentialRate);


#############################################################
#
#	getMeanTipRateStates(....)
#
# returns vector of mean tip states (either net diversification rate or brownian motion
#	rate parameter depending on the 'type' of ephy).
#	use.names=TRUE returns a named vector

getMeanTipRateStates <- function(ephy, use.names = FALSE){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}

	if (use.names == FALSE){
		return(ephy$meanTipLambda - ephy$meanTipMu);
	}
	else{
		ret <- ephy$meanTipLambda - ephy$meanTipMu;
		names(ret) <- ephy$tip.label;
		return(ret);
	}
}

#############################################################
#
#	getStartStopTimes(....)
#
#	adds begin and end times (absolute time) to each edge of 
#	phylogenetic tree

getStartStopTimes <- function(phy){
 	bmax <- max(branching.times(phy));
	bt <- bmax - branching.times(phy);
	begin <- bt[as.character(phy$edge[,1])];
	end <- begin + phy$edge.length;
	phy$begin <- as.numeric(begin);
	phy$end <- as.numeric(end);
	return(phy);
}


#############################################################
#
#	recursivelySetNodeStates(....)
#
# Function to recursively assign states to nodes going from root-to-tip
#	takes args of tree plus vector of node event assignments
#	returns vector.
# phySV is tree with "statevec" component
# node is the focal node
# state is the state (probably "event" in this context)

recursivelySetNodeStates <- function(phySV, node, state) {
	
	phySV$statevec[phySV$edge[,2] == node] <- state;
	if (sum(phySV$edge[,1] == node) > 0){
		# node is internal
		dset <- phySV$edge[,2][phySV$edge[,1] == node];
		phySV <- recursivelySetNodeStates(phySV, dset[1], state);
		phySV <- recursivelySetNodeStates(phySV, dset[2], state);
	}
	
	return(phySV);
}


#############################################################
#
#	getRateThroughTimeMatrix(....)
#	
#	include or exclude options:
#	start.time		=	 start time (in units before present)
#						 if NULL, starts at root
#	end.time		=	 end time 
#						 if NULL, ends at present
#	nslices			=	 number of time cuts
#	Return
#	list with components:
#					$times	= the time vector
#					$lambda = speciation rate matrix
#					$mu 	= extinction rate matrix
# 					$type   = diversification or trait (needs to be extended to trait data)		
# returns object of class bamm-ratematrix
#	 

getRateThroughTimeMatrix <- function(ephy, start.time=NULL, end.time=NULL, nslices=100, node=NULL, nodetype = 'include'){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	if (is.null(node)){
		nodeset <- ephy$edge[,2];
	}else if (!is.null(node) & nodetype == 'include'){
		nodeset <- getDesc(ephy, node)$desc_set;
	}else if (!is.null(node) & nodetype == 'exclude'){
		nodeset <- setdiff( ephy$edge[,2],  getDesc(ephy, node)$desc_set);
	}else{
		stop('error in getRateThroughTimeMatrix\n');
	}
	
	bt <- branching.times(ephy);
		
	maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[,1]))]);

	if (is.null(start.time)){
		start.time <- max(bt) - maxpossible;
	}
	
	if (is.null(end.time)){
		end.time <- max(bt);
	}
 
	tvec <- seq(start.time, end.time, length.out= nslices);
	
	mm <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));
	mumat <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=length(tvec));
	
	for (i in 1:nrow(mm)){
 	
		es <- ephy$eventBranchSegs[[i]];
		events <- ephy$eventData[[i]];
		
		for (k in 1:length(tvec)){
			isGoodTime <- es[,2] <= tvec[k] & es[,3] >= tvec[k];
			if (is.null(node)){ 
				isGoodNode <- rep(TRUE, nrow(es));
			}else{
				isGoodNode <- es[,1] %in% nodeset;	
			}
			estemp <- es[isGoodTime & isGoodNode, ];
			tvv <- tvec[k] - events$time[estemp[,4]];
			rates <- exponentialRate(tvv, events$lam1[estemp[,4]], events$lam2[estemp[,4]]);
			mm[i, k] <- mean(rates);
			mumat[i,k] <- mean(events$mu1[estemp[,4]]);	
		}
		
	}
	
	obj <- list();
	obj$lambda <- mm;
	obj$mu <- mumat;
	obj$times <- tvec;
	
	class(obj) <- 'bamm-ratematrix';
	if(ephy$type=='diversification'){
		obj$type = 'diversification';
	}
	else{
		obj$type = 'traits';	
	}
	return(obj);
}

#############################################################
#
#	getDesc(....)
#
#	 returns a copy of the tree with a component 'desc_set', which
#	is a vector of all desc nodes in downpass sequences order

getDesc <- function(phy, node)
{
	if (is.null(phy$desc_set)){
		phy$desc_set <- node;
	}
	
	if (node > length(phy$tip.label)){
		dset <- phy$edge[,2][phy$edge[,1] == node];
		phy$desc_set <- c(phy$desc_set, dset[1]);
		phy <- getDesc(phy, dset[1]);
		phy$desc_set <- c(phy$desc_set, dset[2]);
		phy <- getDesc(phy, dset[2]);
		
	}
 
 	return(phy);
}

#############################################################
#
#	getSpanningTaxonPair(....)
#
#	returns pair of taxa that span a given taxon set

getSpanningTaxonPair <- function(phy, taxset){
	
	if (! sum(taxset %in% phy$tip.label) > 0){
		cat('Some species in taxset that are not in tree\n');
		taxset <- taxset[taxset %in% phy$tip.label];
	}
	
	dt <- drop.tip(phy, setdiff(phy$tip.label, taxset));
	
	return(c(dt$tip.label[1], dt$tip.label[length(dt$tip.label)]));
}


#############################################################
#
#	getRecursiveSequence(....)
#
#	Private function, called by getEventDataDiversification

getRecursiveSequence <- function(phy){
		
	root.node <- length(phy$tip.label) + 1;
	
	phy$downseq <- root.node
	phy$lastvisit <- numeric(length(unique(phy$edge[,2])));
 
 	phy <- getSequenceForwardTraversal(phy, root.node);

	return(phy);
}

#############################################################
#
#	getRecursiveSequence(....)
#
#	Private function, called by getRecursiveSequence


getSequenceForwardTraversal <- function(phy, node){
	
	if (node <= length(phy$tip.label)){
		#phy$downseq <- c(phy$downseq, node);
	}else{
		dset <- phy$edge[,2][phy$edge[,1] == node];
		phy$downseq <- c(phy$downseq, dset[1]);
		phy <- getSequenceForwardTraversal(phy, dset[1]);
		phy$downseq <- c(phy$downseq, dset[2]);
		phy <- getSequenceForwardTraversal(phy, dset[2]);	
	}
	
	phy$lastvisit[node] <- phy$downseq[length(phy$downseq)];
 
	return(phy);
}

#########

#############################################################
#
#	getCladeRates(....)
#
#	mean clade-specific rates 
#		average of all branch-rates, but weighted by branch length
#		node.type: will compute rates only for clade descended from specified node with 'include'
#					will compute for all branches excluding a given clade, nodetype = 'exclude'
#		

getCladeRates <- function(ephy, weights='branchlengths', node = NULL, nodetype='include', verbose=F){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}	
	
	if (is.null(node)){
		nodeset <- ephy$edge[,2];
	}else if (!is.null(node) & nodetype == 'include'){
		nodeset <- getDesc(ephy, node)$desc_set;
	}else if (!is.null(node) & nodetype == 'exclude'){
		nodeset <- setdiff( ephy$edge[,2],  getDesc(ephy, node)$desc_set);
	}else{
		stop('error in getRateTHroughTimeMatrix\n');
	}
	
	lambda_vector <- numeric(length(ephy$eventBranchSegs));
	mu_vector <- numeric(length(ephy$eventBranchSegs));
 
	
	for (i in 1:length(ephy$eventBranchSegs)){
		if (verbose){
			cat('Processing sample ', i, '\n');
		}
		esegs <- ephy$eventBranchSegs[[i]];
		
		esegs <- esegs[esegs[,1] %in% nodeset, ];
	 
		
		events <- ephy$eventData[[i]];
		events <- events[order(events$index), ];			
		
		# relative start time of each seg, to event:
		relsegmentstart <- esegs[,2] - ephy$eventData[[i]]$time[esegs[,4]];
		relsegmentend <- esegs[,3] - ephy$eventData[[i]]$time[esegs[,4]];
		lam1 <- ephy$eventData[[i]]$lam1[esegs[,4]];
		lam2 <-  ephy$eventData[[i]]$lam2[esegs[,4]];
		mu1 <-  ephy$eventData[[i]]$mu1[esegs[,4]];
		mu2 <-  ephy$eventData[[i]]$mu2[esegs[,4]];
 		
 		seglengths <- esegs[,3] - esegs[,2];	
		wts <- seglengths / sum(seglengths);
		lamseg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, lam1, lam2) / seglengths;
		museg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, mu1, mu2) / seglengths;
	
		lambda_vector[i] <- sum(lamseg * wts);
		mu_vector[i] <- sum(museg  * wts);			
	
	
	}
		
	return(list(lambda = lambda_vector, mu = mu_vector));
}

#############################################################
#
#	timeIntegratedBranchRate(....)
#		computes the integral of rates on a branch segment with respect to time
#		Not the average.
#		Integrates the exponential function p1 * exp(p2 * t)
#		

timeIntegratedBranchRate <- function(t1, t2, p1, p2){
	
	if (p2 == 0){
		return(p1 * (t2 - t1));
	}else{
		(p1/p2)*(exp(p2*t2) - exp(p2*t1));
	}
}
timeIntegratedBranchRate <- Vectorize(timeIntegratedBranchRate);

#########
#############################################################
#
#	getMarginalBranchRateMatrix(....)
#
#	get matrix of marginal rates on each branch for each sample from posterior
# 	
#	BAMM also can directly ouput these marginal rates,
#		so worth having a second function that works directly
#		with the BAMM output file, as it is much faster.
#	getMarginalBranchRateMatrix
getMarginalBranchRateMatrix <- function(ephy, verbose=FALSE){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	lammat <- matrix(0, ncol=length(ephy$eventBranchSegs), nrow=nrow(ephy$edge));
	mumat <- matrix(0, ncol=length(ephy$eventBranchSegs), nrow=nrow(ephy$edge));
	
	for (i in 1:length(ephy$eventBranchSegs)){
		if (verbose){
			cat('Processing sample ', i, '\n');
		}
		esegs <- ephy$eventBranchSegs[[i]];
		events <- ephy$eventData[[i]];
		events <- events[order(events$index), ];			
		
		# relative start time of each seg, to event:
		relsegmentstart <- esegs[,2] - ephy$eventData[[i]]$time[esegs[,4]];
		relsegmentend <- esegs[,3] - ephy$eventData[[i]]$time[esegs[,4]];
		lam1 <- ephy$eventData[[i]]$lam1[esegs[,4]];
		lam2 <-  ephy$eventData[[i]]$lam2[esegs[,4]];
		mu1 <-  ephy$eventData[[i]]$mu1[esegs[,4]];
		mu2 <-  ephy$eventData[[i]]$mu2[esegs[,4]];
 
		lamint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, lam1, lam2);
		muint <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, mu1, mu2);
		seglengths <- esegs[,3] - esegs[,2];	
				
		for (k in 1:nrow(ephy$edge)){
			isRightBranch <- esegs[,1] == ephy$edge[k,2];
			lammat[k, i] <- sum(lamint[isRightBranch]) / sum(seglengths[isRightBranch]);
			mumat[k, i] <- sum(muint[isRightBranch]) / sum(seglengths[isRightBranch]);
			
		}
	
	}
	
	return(list(lambda_branch_matrix = lammat, mu_branch_matrix = mumat));
	
}



#############################################################
#
#	getEventCorrelationMatrix(....)
#
# 	Each entry of this matrix represents the expected 
#	probability that a pair[i, j] of tips will have the same 
#	rate parameters due to BAMM model.
#
# 	Should modify this to allow exponential, spherical,
#		and other possible correlation structures.
#	Need to make a corStruct class that works with this
#		for GLS analyses

getEventCorrelationMatrix <- function(ephy){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	TOL <- 0.0001;
	corMat <- matrix(0, nrow=length(ephy$tip.label), ncol=length(ephy$tip.label));
	for (i in 1:length(ephy$tipStates)){
		dd <- dist(ephy$tipStates[[i]]);
		cmat <- as.matrix(dd);	
		corMat <- corMat + (cmat < TOL);
		
	}
	rownames(corMat) <- ephy$tip.label;
	colnames(corMat) <- ephy$tip.label;
	return(corMat/ length(ephy$numberEvents));
}


#############################################################
#
#	getMeanBranchLengthTree(....)
#
#	Should be capable of computing this from 2 different types of data:
#		1. The bamm-data object
#		2. A list of phylogenetic trees,
#			e.g., as you would get directly as output from BAMM with the 
#			marginal branch length trees.
#		Arg: (for now) 
#			obj	= a list of multiple phylogenetic trees, where each tree
#						has branch length equal to the mean rate on that branch
#						OR a bamm-data object
#			ndr = TRUE or FALSE. if TRUE the mean branch length tree will have the
#				mean net rate on each branch, i.e., lambda - mu. if FALSE the mean
#				rate will simply be lambda.  For trait data ndr will make no difference  

getMeanBranchLengthTree <- function(obj,ndr=TRUE){
	
	if('bamm-data' %in% class(obj))
	{
		v <- list(edge = obj$edge, Nnode = obj$Nnode,tip.label = obj$tip.label, edge.length = obj$edge.length);
		attributes(v) <- list(names = c("edge","Nnode","tip.label","edge.length"),class="phylo",order="cladewise");

		obj <- getMarginalBranchRateMatrix(obj,verbose=FALSE);
		if(ndr)
		{
			el  <- rowMeans(obj$lambda) - rowMeans(obj$mu);
		}
		else
		{
			el  <- rowMeans(obj$lambda);
		}
		v$edge.length <- el;
		obj <- list();
		obj$phy <- v;
		obj$median <- median(el);
		obj$mean <- mean(el);			
	}
	else if(class(obj[[1]]) == 'phylo')
	{
		v <- obj[[1]];
		meanvec <- numeric(length(obj));
 
		el <- v$edge.length;
		meanvec[1] <- mean(v$edge.length);
		for (i in 2:length(obj)){
			el <- el + obj[[i]]$edge.length;
			meanvec[i] <- mean(obj[[i]]$edge.length);
		}
		el <- el / length(obj);
	
		v$edge.length <- el;
		obj <- list();
		obj$phy <- v;
		obj$median <- median(meanvec);
		obj$mean <- mean(meanvec);
	}
	else
	{
		stop("Method not implemented for supplied object class");
	}
	return(obj);		
}

#############################################################
#
#	getSampleDistanceMatrixBAMM <- function(...)
#   For each of k samples from the posterior with branch-specific rates
#		computes pairwise rate differences between sample
#		Thus, for any two rate configurations i and j
#			the i,j and j,i elements of the matrix will be 
#			the squared rate differences
#	
#		phylist is a list of phylogenetic trees with branch-specific rates for 
#			branch lengths (rather than time)

getSampleDistanceMatrixBAMM <- function(phylist, modeltree){
	
	
	# first, scale all edge lengths by time
	for (i in 1:length(phylist)){
		phylist[[i]]$edge.length <- phylist[[i]]$edge.length * modeltree$edge.length;
	}
	
	dmat <- matrix(0, nrow=length(phylist), ncol=length(phylist));
		
	nsamples <- length(phylist);	
	for (i in 1:(nsamples - 1)){
		for (j in i:nsamples){
			dmat[i, j] <- sum((phylist[[i]]$edge.length - phylist[[j]]$edge.length)^2);
		}
	}	
	dmat <- dmat + t(dmat);
	
	return(dmat);
	
}


getSampleCoMat <- function(phylist, modeltree){
	
	comat <- matrix(NA, nrow=length(phylist[[1]]$edge.length), ncol=length(phylist));
	
	# first, scale all edge lengths by time
	for (i in 1:length(phylist)){
		comat[,i] <- phylist[[i]]$edge.length * modeltree$edge.length;
	}	
 
	return(comat);
	
}

#############################################################
#
#	getSpeciesRateThroughTime <- function(...)
#
#	Start time: how many time units before present to include
#	nbreaks: how many time points to compute the rate
#
#	This version computes for individual species
#
getSpeciesRateThroughTime <- function(ephy, start.time, nbreaks=10, ndr=TRUE, species){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
		
	tend <- max(branching.times(ephy))*0.999;
	tstart <- max(branching.times(ephy)) - tend
	tseq <- seq(tstart, tend, length.out=nbreaks);
 
	res <- numeric(nbreaks);
	
	index <- which(ephy$tip.label == species);
	
 
	path <- getPathToRoot(ephy, node=index);
		
	for (k in 1:length(ephy$eventBranchSegs)){

		ed <- ephy$eventData[[k]];
			
		for (z in 1:length(tseq)){
			isGoodNode <- ephy$eventBranchSegs[[k]][,1] %in% path;
			isGoodStart <- ephy$eventBranchSegs[[k]][,2] <= tseq[z];
			isGoodEnd <- ephy$eventBranchSegs[[k]][,3] >= tseq[z];
				
			ev <- ephy$eventBranchSegs[[k]][,4][isGoodNode & isGoodStart & isGoodEnd];
				
			if (ndr){
				lam <- exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
				mu <- ed$mu1[ev];
				res[z] <- res[z] + (lam - mu);	
			}else{
				res[z] <- res[z] + exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
			}
				
				
		}
			
	}
 	res <- res / length(ephy$eventBranchSegs);
 
 	return(res);
	
}


#############################################################
#
#	getBySpeciesRateMatrix <- function(...)
#
#	Start time: how many time units before present to include
#	nbreaks: how many time points to compute the rate
#
#	This version computes for all species
#	
#	node argument will just compute for the subtree descended from "node"

getBySpeciesRateMatrix <- function(ephy, start.time, nbreaks, ndr=TRUE, node){
	
	
	spset <- ephy$tip.label;
	
	mm <- matrix(NA, nrow=length(spset), ncol=nbreaks);
	for (i in 1:nrow(mm)){
		#cat(spset[i], '\n')
		mm[i,] <- getSpeciesRateThroughTime(ephy, start.time, nbreaks, ndr, species=spset[i]);
	}
	
	rownames(mm) <- spset;
	return(mm);
}






#############################################################
#
#	getPathToRoot <- function(...)
#
#	Internal function, gives node path from some node "node" to root
getPathToRoot <- function(phy, node){
	
	root <- length(phy$tip.label) + 1;
	nset <- node;
	while (node != root){
		node <- phy$edge[,1][phy$edge[,2] == node];
		nset <- c(nset, node);
	}
	return(nset);
}


#############################################################
#
#	transparentColor <- function(...)
#
#	Internal function allows for defining named colors with opacity
#	alpha = opacity
#
transparentColor<-function(namedColor,alpha=0.8){
	res<-c(as.vector(col2rgb(namedColor))/255,alpha);
	return(rgb(red=res[1],green=res[2],blue=res[3],alpha=res[4]));
}


#############################################################
#
# plotRateThroughTime <- function(...)
#
# useMedian = boolean, will plot median if TRUE, mean if FALSE.
# intervals if NULL, no intervals will be plotted, otherwise a vector of quantiles must be supplied (these will define shaded polygons)
# ratetype = 'speciation' or 'extinction' or 'netdiv' or 'trait'
# nBins = number of time slices used to generate rates through time
# smooth = boolean whether or not to apply loess smoothing
# smoothParam = loess smoothing parameter, ignored if smooth = F
# opacity = opacity of color for interval polygons
# intervalCol = transparent color for interval polygons
# avgCol = color for mean/median line
# start.time = start time to be fed to getRateThroughTimeMatrix
# end.time = end time to be fed to getRateThroughTimeMatrix
# node = if supplied, the clade descended from this node will be used.
# nodetype = supplied to getRateThroughTimeMatrix
# plot = boolean: if TRUE, a plot will be returned, if FALSE, the data for the plot will be returned. 
#
plotRateThroughTime <- function(ephy, useMedian = F, intervals=seq(from = 0,to = 1,by = 0.01), ratetype = 'speciation', nBins = 100, smooth = F, smoothParam = 0.20, opacity = 0.01, intervalCol='blue', avgCol='red',start.time = NULL, end.time = NULL, node = NULL, nodetype='include', plot = T){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	if (!is.logical(useMedian)){
		stop('ERROR: useMedian must be either TRUE or FALSE.');
	}
	if (class(intervals)!='numeric' & class(intervals)!='NULL'){
		stop("ERROR: intervals must be either 'NULL' or a vector of quantiles.");
	}

	#get rates through binned time
	rmat <- getRateThroughTimeMatrix(ephy, start.time = start.time, end.time = end.time,node = node, nslices = nBins, nodetype=nodetype);

	#set appropriate rates
	if (ratetype != 'speciation' & ratetype != 'extinction' & ratetype != 'netdiv' & ratetype != 'trait'){
		stop("ERROR: ratetype must be 'speciation', 'extinction', 'netdiv' or 'trait'.\n");
	}
	if (ratetype == 'speciation'){
		rate <- rmat$lambda;
		ratelabel <- 'Speciation';
	}
	if (ratetype == 'trait'){
		rate <- rmat$lambda;
		ratelabel <- 'BM rate';
	}

	if (ratetype == 'extinction'){
		rate <- rmat$mu;
		ratelabel <- 'Extinction';
	}
	if (ratetype == 'netdiv'){
		rate <- rmat$lambda - rmat$mu;
		ratelabel <- 'Net diversification';
	}

	#generate coordinates for polygons
	maxTime <- max(rmat$times);
	if (!is.null(intervals)){
		mm <- apply(rate, MARGIN = 2, quantile, intervals);

		poly<-list();
		q1<-1;
		q2<-nrow(mm);
		repeat{
			if (q1 >= q2) {break}
			a<-as.data.frame(cbind(rmat$times,mm[q1,]));
			b<-as.data.frame(cbind(rmat$times,mm[q2,]));
			b<-b[rev(rownames(b)),];
			colnames(a)<-colnames(b)<-c('x','y');
			poly[[q1]]<-rbind(a,b);
			q1<-q1+1;
			q2<-q2-1;
		}
	}

	#Calculate averaged data line
	if (useMedian == F){
		avg <- colMeans(rate);
	}
	if (useMedian == T){
		avg <- unlist(apply(rate,2,median));
	}
	
	#apply loess smoothing to intervals
	if (smooth == T){
		for (i in 1:length(poly)){
			poly[[i]][1:nrow(poly[[i]])/2,2] <- loess(poly[[i]][1:nrow(poly[[i]])/2,2] ~ poly[[i]][1:nrow(poly[[i]])/2,1],span = smoothParam)$fitted;
			poly[[i]][(nrow(poly[[i]])/2):nrow(poly[[i]]),2] <- loess(poly[[i]][(nrow(poly[[i]])/2):nrow(poly[[i]]),2] ~ poly[[i]][(nrow(poly[[i]])/2):nrow(poly[[i]]),1],span = smoothParam)$fitted;
		}
		avg <- loess(avg ~ rmat$time,span = smoothParam)$fitted;
	}

	#begin plotting
	if (plot == T){
		plot.new();
		plot.window(xlim=c(maxTime, 0), ylim=c(0 , max(poly[[1]][,2])));
	
		#plot intervals
		if (!is.null(intervals)){
			for (i in 1:length(poly)){
				polygon(x=maxTime - poly[[i]][,1],y=poly[[i]][,2],col=transparentColor(intervalCol,opacity),border=NA);
			}
		}
		lines(x = maxTime - rmat$time, y = avg, lwd = 3, col = avgCol);

		axis(at=seq(0, maxTime + 0.3*maxTime, by = 5), cex.axis = 1, side = 1);
		axis(at=seq(-0.2, max(rate) + 0.2*max(rate), by=0.1), las=1, cex.axis = 1, side = 2);
		mtext(side = 1, text = 'Time since present', line = 3, cex = 1.1);
		mtext(side = 2, text = ratelabel, line = 3, cex = 1.1);
	}
	if (plot == F){
		return(list(poly = poly,avg = avg,times = rmat$time))
	}
}


##################################################################
#
#	getSpeciesRateThroughTimeReturnAll <- function(...)
#
#	Primarily intended as an internal function to support plotSpeciesRatesThroughTime()
#	Identical to getSpeciesRateThroughTime() except that this function returns all values from posterior, rather than averaging.
#
#
getSpeciesRateThroughTimeReturnAll <- function(ephy, start.time, nbreaks, ndr, species){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
		
	tend <- max(branching.times(ephy))*0.999;
	tstart <- max(branching.times(ephy)) - tend;
	tseq <- seq(tstart, tend, length.out=nbreaks);
 
	res <- numeric(nbreaks);
	
	index <- which(ephy$tip.label == species);
	
	path <- getPathToRoot(ephy, node=index);
	
	resMat <- matrix(NA, nrow=length(ephy$eventBranchSegs), ncol=nbreaks);
	for (k in 1:length(ephy$eventBranchSegs)){

		ed <- ephy$eventData[[k]];
		
		resVec <- vector(mode='numeric', length=nbreaks);
		for (z in 1:length(tseq)){
			isGoodNode <- ephy$eventBranchSegs[[k]][,1] %in% path;
			isGoodStart <- ephy$eventBranchSegs[[k]][,2] <= tseq[z];
			isGoodEnd <- ephy$eventBranchSegs[[k]][,3] >= tseq[z];
				
			ev <- ephy$eventBranchSegs[[k]][,4][isGoodNode & isGoodStart & isGoodEnd];
				
			if (ndr){
				lam <- exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
				mu <- ed$mu1[ev];
				resVec[z] <- lam - mu;	
			}else{
				resVec[z] <- exponentialRate(tseq[z]-ed$time[ev], ed$lam1[ev], ed$lam2[ev]);
			}	
		}		
	resMat[k,] <- resVec;
	}
	return(resMat);
}



#########################################################
#
#	plotSpeciesRatesThroughTime <- function(...)
#
#	Function to plot several species-specific rate trajectories together
#
#	start.time = time before present (defaults to full depth of tree)
#	useMedian = boolean, will plot median if TRUE, mean if FALSE
#	nbreaks = number of time slices to use
#	ratetype = 'speciation' or 'netdiv' or 'trait'
#	species = vector of species names to include
#	lowerCI = lowest quantile to consider when plotting confidence intervals
#	upperCI = greatest quantile to consider when plotting confidence intervals
#		if either upperCI or lowerCI are NULL, then confidence intervals will be ignored.
#	smooth = boolean, apply loess smoothing to curves
#	spCol = vector of colors that will be matched to the species vector
#	opacity = opacity level for plotting colored confidence intervals


plotSpeciesRatesThroughTime <- function(ephy, start.time=max(branching.times(ephy)), useMedian=F, nbreaks=10, ratetype='speciation', species, lowerCI=0.05, upperCI=0.95, smooth=T, spCol, opacity=0.01,smoothParam=0.20){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	if (!all(species %in% ephy$tip.label)){
		stop('Not all species names match those in ephy.');
	}
	
	if (ratetype=='speciation' | ratetype=='trait'){
		ndr <- F;
		if (ratetype == 'speciation'){ratelabel <- 'Speciation'}
		if (ratetype == 'trait'){ratelabel <- 'BM rate'}
	}
	if (ratetype == 'netdiv'){
		ndr <- T;
		ratelabel <- 'Net diversification';
	}
		
	bySp<-lapply(species,function(x) getSpeciesRateThroughTimeReturnAll(ephy, start.time = start.time, nbreaks = nbreaks, ndr = ndr, species = x));
	names(bySp) <- species;
	
	#calculate average values
	if (useMedian == T){
		avgList<-lapply(bySp,function(x) apply(x,2,function(x) median(x)));
	}
	if (useMedian == F){
		avgList<-lapply(bySp,function(x) colMeans(x));
	}
	
	tend <- max(branching.times(ephy))*0.999;
	tstart <- max(branching.times(ephy)) - tend;
	tseq <- seq(tstart, tend, length.out = nbreaks);
	
	#calculate confidence intervals
	if (is.numeric(lowerCI) & is.numeric(upperCI)){
		intervals <- seq(from = lowerCI,to = upperCI,by = 0.01);
		
		polySp<-list();
		for (i in 1:length(bySp)){
			mm <- apply(bySp[[i]], 2, quantile, intervals);
			poly <- list();
			q1 <- 1;
			q2 <- nrow(mm);
			repeat{
				if (q1 >= q2) {break}
				a <- as.data.frame(cbind(tseq,mm[q1,]));
				b <- as.data.frame(cbind(tseq,mm[q2,]));
				b <- b[rev(rownames(b)),];
				colnames(a) <- colnames(b) <- c('x','y');
				poly[[q1]] <- rbind(a,b);
				q1 <- q1+1;
				q2 <- q2-1;
			}
			polySp[[i]] <- poly;
		}
	}
	
	if (smooth == T){
		if (nbreaks < 30){
			cat('Too few breaks. Non-smoothed results returned.\n');
		}
		if (nbreaks >= 30){
			for (i in 1:length(polySp)){
				for (j in 1:length(polySp[[i]])){
					polySp[[i]][[j]][1:nrow(polySp[[i]][[j]])/2,2] <- loess(polySp[[i]][[j]][1:nrow(polySp[[i]][[j]])/2,2] ~ polySp[[i]][[j]][1:nrow(polySp[[i]][[j]])/2,1],span = smoothParam)$fitted;
					polySp[[i]][[j]][(nrow(polySp[[i]][[j]])/2):nrow(polySp[[i]][[j]]),2] <- loess(polySp[[i]][[j]][(nrow(polySp[[i]][[j]])/2):nrow(polySp[[i]][[j]]),2] ~ polySp[[i]][[j]][(nrow(polySp[[i]][[j]])/2):nrow(polySp[[i]][[j]]),1],span = smoothParam)$fitted;
				}
			}
			for (i in 1:length(avgList)){
				avgList[[i]] <- loess(avgList[[i]] ~ tseq,span = smoothParam)$fitted;
			}
		}
	}
	
	maxRate <- max(unlist(lapply(polySp,function(x) lapply(x,function(x) x[,2]))));
	plot.new();
	plot.window(xlim=c(max(branching.times(ephy)), 0), ylim=c(0 , maxRate));

	if (is.numeric(lowerCI) & is.numeric(upperCI)){
		for (i in 1:length(polySp)){
			for (j in 1:length(polySp[[i]])){
				polygon(x = tend - polySp[[i]][[j]][,1],y = polySp[[i]][[j]][,2], col = transparentColor(spCol[i],alpha = opacity), border = NA);
			}
			lines(x = tend - tseq,y = avgList[[i]], lwd=2, col = spCol[i]);
		}
	}
	if (!is.numeric(lowerCI) | !is.numeric(upperCI)){
		for (i in 1:length(avgList)){
			lines(x = tend - tseq,y = avgList[[i]],lwd = 2,col = spCol[i]);
		}
	}
	
	axis(at = seq(0, + 1.3*max(branching.times(ephy)), by = 5), cex.axis = 1, side = 1);
	axis(at = seq(-0.2, maxRate + 0.2*maxRate, by=0.2), las=1, cex.axis = 1, side = 2);
	mtext(side = 1, text = 'Time since present', line = 3, cex = 1.1);
	mtext(side = 2, text = ratelabel, line = 3, cex = 1.1);
}


#############################################################
#
#	marginalShiftProbsTree(....)
#
#	Args: ephy	=	object of class 'bamm-data'
#	
#	Returns:		a phylogenetic tree, but where each 
#	             	branch length (edge length) is equal to the
#					marginal probability of shift occuring 
#					on that particular branch. 
#							

marginalShiftProbsTree <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}
	
	shiftvec <- numeric(length(ephy$edge.length));
 	rootnode <- length(ephy$tip.label) + 1;
 
 
	for (i in 1:length(ephy$eventData)){
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		shiftvec[hasShift] <- shiftvec[hasShift] + rep(1, sum(hasShift));
	}
	
	shiftvec <- shiftvec / length(ephy$eventData);	
	
	newphy <- as.phylo.bammdata(ephy);
	newphy$edge.length <- shiftvec;
	return(newphy);
}


#############################################################
#
#	cumulativeShiftProbsTree(....)
#
#	Args: ephy	=	object of class 'bamm-data'
#	
#	Returns:		a phylogenetic tree, but where each 
#	             	branch length (edge length) is equal to the
#					cumulative probability of a shift somewhere 
#					between the focal branch and the root of the 
#					tree. The branch length itself does not tell 
#					you where the shifts occur, but they tell 
#					you which clades/lineages have diversification
#					dynamics that are decoupled from the root of the tree 
#							

cumulativeShiftProbsTree <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}	
			
	shiftvec <- numeric(length(ephy$edge.length));
 	rootnode <- length(ephy$tip.label) + 1;
 
 
	for (i in 1:length(ephy$eventData)){
		
		snodes <- unique(ephy$eventBranchSegs[[i]][,1][ephy$eventBranchSegs[[i]][,4] != 1]);
		hasShift <- ephy$edge[,2] %in% snodes;
		shiftvec[hasShift] <- shiftvec[hasShift] + rep(1, sum(hasShift));
	}	
	
	shiftvec <- shiftvec / length(ephy$eventData);		
	newphy <- as.phylo.bammdata(ephy);
	newphy$edge.length <- shiftvec;
	return(newphy);	
	
	
}

#############################################################
#
#	maximumShiftCredibilityTree(....)
#
#	Args: ephy	=	object of class 'bamm-data'
#	
#
#	maximize		=	'sum', = sum of branch probabilities for each tree
#						'product', = product of branch probabilities for each tree
#	
#	Returns: 		- bestconfigs: a list of length equal the number of 
#						unique shift configurations in the maximum shift
#						credibility set. Each element is a vector of sample
#						indices from the BAMM-data object with identical
#						shift configurations.
#					  
#					- A vector of optimality scores for all other samples 
#						in posterior from the bamm-data object.
#
#					- sampleindex: a representative index for samples from each 
#						set of unique shift configurations. The length of this vector
#						is equal to the length of the bestconfigs list. If this vector was
#						sampleindex = c(2, 20, 50), this would mean that there are 3 distinct
#						sets of shift configurations with equal credibility under the optimality 
#						criterion. More commonly, a single shift configuration will be dominant, and 
#						although the length of bestconfigs[[1]] may be greater than 1, the sampleindex
#						vector will contain a single representative event from that set.
#
#	See example file.
#   This is analogous to the maximum clade credibility tree from a 
#		Bayesian phylogenetic analysis.

maximumShiftCredibilityTree <- function(ephy, maximize = 'product'){

	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}			
	
	probvec <- numeric(length(ephy$eventData));
	
	mtree <- marginalShiftProbsTree(ephy);
	px <- mtree$edge.length;
	

	for (i in 1:length(ephy$eventData)){
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		branchprobs <- (hasShift)*px  + (!hasShift)*(1 - px) ;
		if (maximize == 'product'){
			probvec[i] <- sum(log(branchprobs));
		}else if (maximize == 'sum'){
			probvec[i] <- sum(branchprobs);
		}else{
			stop("Unsupported optimize criterion in maximumShiftCredibilityTree");
		}
	}
	
	best <- which(probvec == max(probvec));
	
	# Now test for multiple trees with same log-prob:
	bestconfigs <- list();
		
	index <- 0;	
	while (length(best) > 0){
		index <- index + 1;	
		lv <- logical(length = length(best));
		for (i in 1:length(best)){
			lv[i] <- areEventConfigurationsIdentical(ephy, best[1], best[i]);
		}
		bestconfigs[[index]] <- best[lv];
		best <- best[!lv];
	}
	
	
	sampleindex <- numeric(length(bestconfigs));
	for (i in 1:length(bestconfigs)){
		sampleindex[i] <- bestconfigs[[i]][1];
	}
	
	obj <- list();
	obj$bestconfigs <- bestconfigs;
	obj$scores <- probvec;
	obj$optimalityType = maximize;
 	obj$sampleindex <- sampleindex;
	return(obj);
}

areEventConfigurationsIdentical <- function(ephy, index1, index2){
	
	nodeset <- c(ephy$eventData[[index1]]$node, ephy$eventData[[index2]]$node);
	diffs <- sum(table(nodeset) == 1);
	return(diffs == 0);	
}

#############################################################
#
#	getShiftNodesFromIndex (....)
#
#	Args: ephy	=	object of class 'bamm-data'
#		  index =   the index of the sample you wish to view, e.g.,
#					 if index = 5, this will give you the nodes subtending
#                    all branches with rate shifts for the 5th sample
#					 from the posterior in your BAMM-data object.								 
#
#	
#	Returns: 		- a vector of the nodes where shifts occurred, excluding the root.
#					Note: if NO shifts occured, this will return a 
#							numeric object of length zero
# 

getShiftNodesFromIndex <- function(ephy, index){
	
	if (index > length(ephy$eventData)){
		cat("Error. Attempting to access non-existent element from BAMM-data object\n");
		cat("You have << ", length(ephy$eventData), " >>> samples in your BAMM-data object\n");
		cat("Value of index must be no greater than the number of samples\n\n");
		stop();
	}
	root <- length(ephy$tip.label) + 1;
	nodes <- ephy$eventData[[index]]$node;
	nodes <- nodes[nodes != root];

	return(nodes);
}


as.phylo.bammdata <- function(ephy){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}		
	
	newphylo <- list();
	newphylo$edge <- ephy$edge;
	newphylo$Nnode <- ephy$Nnode;
	newphylo$tip.label <- ephy$tip.label;
	newphylo$edge.length <- ephy$edge.length;
	class(newphylo) <- 'phylo';
	return(newphylo);
}


#	type		= the type of analysis to autotune (speciationextinction or trait)
#					should probably correspond to default bamm modeltype options
#	burn		=	fraction of samples from the tuning output file to discard 
#						(recommend 0.1)
#	outfile		= 	file where summary will be stored.
#	
#	opt_prob	= theoretical optimum mcmc acceptance rate
#				0.234 is probably too low, but 
#				anything between 0.23 and 0.4 is probably fine.
#				(0.234 being optimum for infinite-dimensional target space)
#	
#
#
# Parameters by index, for diversification:
#		0 = change number of events (not tuned)
#		1 = event position scale (event location)
#		2 = updateEventRateScale
#		3 = updateLambdaInitScale
#		4 = updateLambdaShiftScale
#		5 = updateMuInitScale
#	

autotune_BAMM_MCMC <- function(fname, type = 'speciationextinction', burn=0.1, opt_prob = 0.234, outfile = 'Autotune_results.txt', min.set = 0.01){
	
	xx <- read.csv(fname, header=F);
	mn <- floor(0.1*nrow(xx));
	xx <- xx[mn:nrow(xx), ];
	
 	if (type == 'speciationextinction'){
 		
 		parnames = c('updateEventLocationScale', 'updateEventRateScale', 
 				'updateLambdaInitScale', 'updateLambdaShiftScale', 'updateMuInitScale');
 		
 		opt_vals <- numeric(length(parnames));
 		comments <- character(length(parnames));
 		
 		for (i in 1:5){
 			xtemp <- xx[xx[,1] == i, ];
 			fit_temp <- glm(xtemp[,3] ~ xtemp[,2], family = 'binomial');
 			opt_vals[i] <- (logit(opt_prob) - fit_temp$coefficients[1]) / fit_temp$coefficients[2];
 			if (opt_vals[i] < 0){
 				comments[i] <- ' Warning: this parameter not tuned. Setting to minimum value';
 				opt_vals[i] <- min.set;
 			}
 		
 		}		
 		
 		setup_outfile(outfile, parnames, opt_vals, comments);	
 		
 	}else if (type == 'trait'){
 		
 		cat("trait auto-tune not yet supported\n");
 	}
	
	
	
}


setup_outfile <- function(fname, parnames, parvalues, comments){
	cat('#\n#\tFile generated ', as.character(Sys.time()), '\n', file = fname, append=F);
	cat('#\t\t\tby BAMMtools::autotune_BAMM_MCMC\n', file = fname, append=T);	
	cat('#\n#\tThis file contains autotuned parameters\n', file = fname, append=T);
	cat('#\tfor a BAMM analysis\n#\n', file = fname, append=T);
	cat('#\tCopy the block of parameter below\n', file = fname, append=T);	
	cat('#\tinto the MCMC_OPERATORS section of your controlfile\n', file = fname, append=T);		
	cat('#\tEnsure that the default MCMC_SCALING_OPERATORS section is deleted\n', file = fname, append=T);	
	cat('#\tstart MCMC_SCALING_OPERATORS block <<<<\n\n', file = fname, append=T);	
	for (i in 1:length(parnames)){
		cat(parnames[i], ' = ', parvalues[i], sep='', file = fname, append=T);
		if (comments[i] != ''){
			cat('   #', comments[i], sep='', file = fname, append=T);				
		}
		cat('\n', file = fname, append=T);		
	}
	cat('\n#\t>>>>>end MCMC_SCALING_OPERATORS block', file = fname, append=T);	
	cat('\n#\tNote: if you received a warning message after a parameter,', file = fname, append=T);	
	cat('\n#\tthis just means that the auto-tuning could not identify a best target value', file = fname, append=T);	
	cat('\n#\tfor this scaling operator. You do not need to worry about this.', file = fname, append=T);	
	cat('\n#\tHowever, if MCMC is performing poorly, you may consider manuually', file = fname, append=T);	
	cat('\n#\tIncreasing or decreasing the value of the parameter', file = fname, append=T);	
}

inv.logit <- function (x, min = 0, max = 1) 
{
    p <- exp(x)/(1 + exp(x))
    p <- ifelse(is.na(p) & !is.na(x), 1, p)
    p * (max - min) + min
}

logit <- function (x, min = 0, max = 1) 
{
    p <- (x - min)/(max - min)
    log(p/(1 - p))
}







#############################################################
#	computeBayesFactors
#
#
#   postfilename		=   MCMC output file from regular BAMM run
#				 			e.g., with sampleFromPriorOnly = 0
#			
#	priorfilename		=	MCMC output file from running BAMM with 
#                           sampleFromPriorOnly = 1
#
#
#	burnin				=	How many samples from posterior to discard
#							Uses fraction (e.g, 0.25 = 25% discarded)
#							Will also discard this same fraction from the prior.
#			
#	modelset			=	Integer set corresponding to models for which
#							you wish to compute pairwise Bayes Factors
#							e.g., 0:2 will compute all pairwise BF between models 
#							with 0 to 2 process 
#							(0 is a model with zero non-root processes)
#	
#							Will only compute Bayes Factors for the set of models
#							0:K that includes 99.5% of the sampled models. 
#
#   Returns:  matrix w pairwise Bayes Factors
#	By convention, the model with the higher index is the numerator for the calculation
#	e.g., M2 / M1 or M1 / M0, but never M0 / M1.
#
#	By default, odds ratios are computed as 
#			(prior_odds_M2  + 1   ) / (prior_odds_M1 + 1)
#		     where the 1 is added to both numerator and denominator 
#			 to avoid divide by zero erros 
#	
 
	
computeBayesFactors <- function(postfilename, priorfilename, burnin = 0.1, modelset = 0:5){

	if (length(modelset) < 2){
		stop('\nInvalid modelset argument. This must be a vector of length > 1');
	}
	
	
	post <- read.csv(postfilename, header=T);
	prior <- read.csv(priorfilename, header=T);
	
	post <- post[floor(burnin*nrow(post)):nrow(post), ];
	prior <- prior[floor(burnin*nrow(prior)):nrow(prior), ];

	
	tpost <- table(post$numevents);
	tprior <- table(prior$numevents);

	fprobs <- cumsum(tpost) / sum(tpost);
 	max_model <- NA;
 	if (length(fprobs) == 1){
 		max_model <- as.numeric(names(fprobs));
 	}else{
 		fprobs <- fprobs[fprobs < 0.995]; 		
 	 	max_model <- as.numeric(names(fprobs[length(fprobs)]));
	}
	
	if (max_model < max(modelset)){
		cat('*****************************************\n');
		cat('You have selected to compute Bayes Factors for models');
		cat('\n that were sampled very infrequently and for which\n');
		cat(' the Bayes Factors are likely to be (wildly ) inaccurate.\n');
		cat(' Consequently, the maximum rank of the models considered\n');
		cat(' will be constrained to <<< ', max_model, ' >>>\n');
		cat('*****************************************\n\n');

	}
	
	modelset <- modelset[modelset <= max_model];
	mset <- as.character(modelset);

	
	postf <- numeric(length(modelset));
	names(postf) <- mset;
	inboth <- intersect(mset, names(tpost));
	postf[inboth] <- tpost[inboth];
	
	priorf <- numeric(length(modelset));
	names(priorf) <- mset;
	inboth <- intersect(mset, names(tprior));
	priorf[inboth] <- tprior[inboth];
	
	mm <- matrix(NA, nrow=length(mset), ncol=length(mset));
	rownames(mm) <- mset;
	colnames(mm) <- mset;
	
	if (length(modelset) < 2){
		stop('\nError. Invalid model choice - is rank of specified model too high?\n');
	}
	
	for (i in 1:(length(modelset) - 1)){
		for (j in (i+1):length(modelset)){
					
			prior_odds <- (priorf[j] + 1) / (priorf[i] + 1);
			post_odds <- (postf[j] + 1) / (postf[i] + 1);
			
			mm[i , j] <- post_odds / prior_odds;
			
		}	
		
	}
	
	return(mm);
	
}











