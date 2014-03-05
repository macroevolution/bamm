
# Feb 28 2014
credibleShiftSet <- function(ephy, threshold, set.limit=0.95){
	
	dsc <- distinctShiftConfigurations(ephy, threshold);
	cfreq <- cumsum(dsc$frequency);
	cut <- min(which(cfreq >= set.limit));
	nodeset <- NULL;

	if (class(threshold) == 'branchprior'){
		dsc$marg.probs <- dsc$marg.probs[as.character(threshold$edge[,2])];	
		nodeset <- dsc$marg.probs[dsc$marg.probs >= threshold$edge.length];
	}else if (class(threshold) == 'numeric'){
		nodeset <- dsc$marg.probs[dsc$marg.probs >= threshold];		
	}else{
		stop("arg threshold is of the wrong class\n");
	}

 	
 	shiftnodes <- dsc$shifts[1:cut];
	indices <- dsc$samplesets[1:cut];
	frequency <- dsc$frequency[1:cut];
	cumulative <- cumsum(dsc$frequency)[1:cut];
	
	root <- length(ephy$tip.label) + 1;
	
	eventVectors <- vector("list", length(indices));
	eventData <- vector("list", length(indices));
	tipStates <- vector("list", length(indices));
	eventBranchSegs <- vector("list", length(indices));
	numberEvents <- vector("list", length(indices));
	tipLambda <- vector("list", length(indices));
	tipMu <- vector("list", length(indices));
		
 	stoptime <- max(ephy$end);	
	
	for (i in 1:length(indices)) {
		
		subd <- ephy$eventData[unlist(indices[i])];
		goodnodeset <- c(root, unlist(shiftnodes[[i]]));
		
		tmp <- subd[[1]];
		#tmp <- tmp[tmp$node %in% goodnodeset, ];
		#tmp <- tmp[order(tmp$time), ];
		
		if (length(subd) > 1){
 			for (k in 2:length(subd)){
 				tmp <- rbind(tmp, subd[[k]]);
 			}
	 
		}
		tmp <- tmp[tmp$node %in% goodnodeset, ];
		
		zz <- numeric(length(goodnodeset));
		dftemp <- data.frame(node=zz, time=zz, lam1=zz, lam2=zz, mu1=zz, mu2=zz, index=zz);
		dftemp$node <- goodnodeset;
		
		for (k in 1:nrow(dftemp)){
			subs <- tmp[tmp$node == dftemp$node[k], ];
			dftemp$time[k] <- mean(subs$time);
			dftemp$lam1[k] <- mean(subs$lam1);
			dftemp$lam2[k] <- mean(subs$lam2);
			dftemp$mu1[k] <- mean(subs$mu1);
			dftemp$mu2[k] <- mean(subs$mu2);
		}
		dftemp <- dftemp[order(dftemp$time ), ];
		dftemp$index <- 1:nrow(dftemp);
		
		rownames(dftemp) <- NULL;
		
		statevec <- rep(1, nrow(ephy$edge));

		if (nrow(dftemp) > 1) {
			for (k in 2:nrow(dftemp)) {
				s1 <- which(ephy$downseq == dftemp$node[k]);
				s2 <- which(ephy$downseq == ephy$lastvisit[dftemp$node[k]]);
				descSet <- ephy$downseq[s1:s2];
				isDescendantNode <- ephy$edge[,2] %in% descSet;				
				statevec[isDescendantNode] <- k;
			}				
		}

 		tmpEventSegMat <- matrix(0, nrow=(max(ephy$edge) + nrow(dftemp) - 2), ncol=4);
 		
		non.root <- c(1:length(ephy$tip.label), (length(ephy$tip.label)+2):max(ephy$edge));
		pos <- 1;	
		
		is_noEventBranch = !(ephy$edge[,2] %in% dftemp$node);
		
		tmpEventSegMat[1:sum(is_noEventBranch), 1] <- ephy$edge[,2][is_noEventBranch];
		tmpEventSegMat[1:sum(is_noEventBranch),2] <- ephy$begin[is_noEventBranch];
 		tmpEventSegMat[1:sum(is_noEventBranch),3] <- ephy$end[is_noEventBranch];
 		tmpEventSegMat[1:sum(is_noEventBranch),4] <- statevec[is_noEventBranch];		
 		
		eventnodeset <- intersect(non.root, dftemp$node);
		pos <- 1 + sum(is_noEventBranch);
		for (k in eventnodeset) {
			events.on.branch <- dftemp[dftemp$node == k, ];
			events.on.branch <- events.on.branch[order(events.on.branch$time), ];
				
			fBranch <- ephy$edge[,2] == k;
 			start.times <- c(ephy$begin[fBranch], events.on.branch$time);
			stop.times <- c(events.on.branch$time, ephy$end[fBranch]);
			parent <- ephy$edge[,1][ephy$edge[,2] == k];
			if (parent == (length(ephy$tip.label) + 1)) {
				# Parent is root:
				proc.set <- c(1, events.on.branch$index);	
			} else {
				proc.set <- c(statevec[ephy$edge[,2] == parent], events.on.branch$index);			
			}
				
 			zzindex = pos:(pos+nrow(events.on.branch));	
				
			tmpEventSegMat[zzindex, 1] <- rep(k, length(zzindex));
			tmpEventSegMat[zzindex, 2] <- start.times;
			tmpEventSegMat[zzindex, 3] <- stop.times;
			tmpEventSegMat[zzindex, 4] <- proc.set;		
			pos <- pos + 1 + nrow(events.on.branch);
		}
		
 		tmpEventSegMat <- tmpEventSegMat[order(tmpEventSegMat[,1]), ];
 	
 		eventBranchSegs[[i]] <- tmpEventSegMat;

		tipstates <- numeric(length(ephy$tip.label));
		tipstates <- statevec[ephy$edge[,2] <= ephy$Nnode+1];
		tipstates <- tipstates[order(ephy$edge[ephy$edge[,2] <= ephy$Nnode+1,2])];
		
 		### Compute tip rates:
		
		#tiplam <- dftemp$lam1[tipstates] * exp(dftemp$lam2[tipstates] * (stoptime - dftemp$time[tipstates]));
		tiplam <- exponentialRate(stoptime - dftemp$time[tipstates], dftemp$lam1[tipstates], dftemp$lam2[tipstates]);
		tipmu <- dftemp$mu1[tipstates];

		### List assignments and metadata across all events:
		eventData[[i]] <- dftemp;	
		eventVectors[[i]] <- statevec;
		numberEvents[i] <- nrow(dftemp);
		tipStates[[i]] <- tipstates;
		
		tipLambda[[i]] <- tiplam;
		tipMu[[i]] <- tipmu;	
 	}
 	
 	obj <- as.phylo.bammdata(ephy);
 	
	obj$begin <- ephy$begin;
	obj$end <- ephy$end; 	
 	
 	obj$downseq <- ephy$downseq;
 	obj$lastvisit <- ephy$lastvisit;
 	
 	obj$numberEvents <- numberEvents;
	obj$eventData <- eventData;
	obj$eventVectors <- eventVectors;
	obj$tipStates <- tipStates;
	obj$tipLambda <- tipLambda;
	obj$meanTipLambda <- rep(NA, length(ephy$meanTipLambda));
	

	obj$eventBranchSegs <- eventBranchSegs; 	
	obj$tipMu <- tipMu;
	obj$meanTipMu <- rep(NA, length(ephy$meanTipMu));
	
	if (ephy$type == 'diversification') {
		obj$type <- 'diversification';
	} 
	else if (ephy$type == 'trait') {
		obj$type <- 'trait';	
	}else{
		stop('Problem: ephy of wrong type');
	}
	
	obj$marg.probs <- dsc$marg.probs;
 	
 	obj$shiftnodes <- shiftnodes;
 	obj$indices <- indices;
 	obj$frequency <- frequency;
 	obj$cumulative <- cumulative;
 	obj$threshold <- threshold;
 	obj$set.limit <- set.limit;
 	obj$number.distinct <- length(indices);
 	
	class(obj) <- 'credibleshiftset';
	return(obj);	
	
}






