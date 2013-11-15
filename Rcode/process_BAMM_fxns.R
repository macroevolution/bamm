
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
getEventData <- function(phy, eventfilename, burnin=0, nsamples = NULL, verbose=F, assign.type = 'new_way', type = 'diversification', header=TRUE){
	
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
 	x <- read.csv(eventfilename, header=header, stringsAsFactors=F);
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
		dftemp <- data.frame(node=nodeVec, time=tm, lam1=lam1, lam2=lam2, mu1=mu1, mu2=mu2, stringsAsFactors=F);
		
		
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
#	getMeanNetDivRate(....)
#
# returns matrix with mean net diverisification rates at tip of tree
#	use.names=T returns a named vector

getMeanNetDivRate <- function(ephy, use.names = FALSE){
	
	if (!'bamm-data' %in% class(ephy)){
		stop("Object ephy must be of class bamm-data\n");
	}

	if (use.names == F){
		return(ephy$meanTipLambda - ephy$meanTipMu);
	}
	if (use.names == T){
		tmp <- ephy$meanTipLambda - ephy$meanTipMu;
		names(tmp) <- ephy$tip.label;
		return(tmp);
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
		stop('error in getRateTHroughTimeMatrix\n');
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
#			phylist	= a list of multiple phylogenetic trees, where each tre
#						has branch length equal to the mean rate on that branch
#		

getMeanBranchLengthTree <- function(obj,ndr=TRUE){
	
	if('bamm-data' %in% class(obj))
	{
		v <- list(edge = obj$edge, Nnode = obj$Nnode,tip.label = obj$tip.label, edge.length = obj$edge.length);
		attributes(v) <- list(names = c("edge","Nnode","tip.label","edge.length"),class="phylo",order="cladewise");

		obj <- getMarginalBranchRateMatrix(obj);
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
	tstart <- start.time #tend - start.time;
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













