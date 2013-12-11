getEventData = function(phy, eventdata, burnin=0, nsamples = NULL, verbose=FALSE, type = 'diversification', header=TRUE)
{	
	if (type != 'diversification' & type != 'trait')
	{
		stop("Invalid 'type' specification. Should be 'diversification' or 'trait'");
	}
	
	if (any(is.null(c(phy$begin, phy$end))))
	{
		phy = getStartStopTimes(phy);
	}
		
	bt = branching.times(phy);
	
	eventVectors = list();
	eventData = list();
	tipStates = list();
	eventBranchSegs = list();
	
	tipLambda = list();
 
	if (class(eventdata) == 'data.frame')
	{
		cat("Processing event data from data.frame\n");
		uniquegens = sort(unique(eventdata[,1]));
	}
	else if (class(eventdata) == 'character')
	{
		cat("Reading event datafile: ", eventdata, "\n\t\t...........");
		eventdata = read.csv(eventdata, header=header, stringsAsFactors=FALSE);
 		uniquegens = sort(unique(eventdata[,1]));
 		cat("\nRead a total of ", length(uniquegens), " samples from posterior\n");				
	}
	else
	{
		err.string = c('eventdata arg invalid\n\nType is ', class(eventdata), '\n', sep='');
		stop(err.string);
	}
 	
 	samplestart = uniquegens[floor(burnin*length(uniquegens))];
 	if(!length(samplestart))
 	{
 		samplestart = 0;
 	}
 	uniquegens = uniquegens[uniquegens >= samplestart];
 
 	if (is.null(nsamples))
 	{
 		nsamples = length(uniquegens);
 	}
 	else if (nsamples > length(uniquegens))
 	{
 		nsamples = length(uniquegens);
 	}
	
	goodsamples = uniquegens[seq.int(1, length(uniquegens), length.out=nsamples)];
 	 
 	cat('\nDiscarded as burnin: GENERATIONS < ', goodsamples[1]);
 	cat("\nAnalyzing ", length(goodsamples), " samples from posterior\n");

 	numberEvents = length(goodsamples); # vector to hold number of events per sample
 
 	cat('\nSetting recursive sequence on tree...\n');
 	phy = getRecursiveSequence(phy);
 	cat('\nDone with recursive sequence\n\n');
 
	######### Get ancestors for each unique pair of taxa
	if (verbose)
	{
		cat("Start preprocessing unique MRCA pairs....\n");
	}	
		
	x2 = eventdata[eventdata$generation %in% goodsamples, ];
		
	uniquePairSet = matrix(NA, nrow=nrow(x2), ncol=2);	
	uniquePairNode = numeric(nrow(x2));
	
	uniquePairSet[,1] = as.integer(match(x2$leftchild, phy$tip.label));
	uniquePairSet[,2] = as.integer(match(x2$rightchild, phy$tip.label, nomatch = 0L));
	uniquePairNode = getmrca(phy, uniquePairSet[,1],uniquePairSet[,2]);
	
	if (verbose)
	{
		cat("Done preprocessing unique MRCA pairs....\n");
	}	

	####### Done with risky sstuff
	
	meanTipMu = numeric(length(phy$tip.label));
	
 	meanTipLambda = numeric(length(phy$tip.label)); 
 
 	for (i in 1:length(goodsamples))
 	{	
  		tmpEvents = x2[x2[,1] == goodsamples[i], ];
		
		if (verbose) cat('Processing event: ', i, '\n');		
 
 		tm = tmpEvents[,4]; # abs time of event
 		lam1 = tmpEvents[,5]; # lambda parameter 1
 		lam2 = tmpEvents[,6]; # lambda parameter 2
 		if(type == 'diversification')
 		{	
			mu1 = tmpEvents[, 7]; # mu parameter 1
 			mu2 = tmpEvents[, 8]; #mu parameter 2 
		}
		else 
		{ #for bamm trait data we set the mu columns to zero because those params don't exist	
			mu1 = rep(0, nrow(tmpEvents)); 
 			mu2 = rep(0, nrow(tmpEvents)); 
		}
		tipMu = list();	
		
 		# Get subtending node for each event:
 		nodeVec = uniquePairNode[x2[,1] == goodsamples[i]];
 		
 		if (sum(nodeVec == 0)  > 0)
 		{
			stop('Failed to assign event to node\n');
		}
		
		# make a dataframe:
		dftemp = data.frame(node=nodeVec, time=tm, lam1=lam1, lam2=lam2, mu1=mu1, mu2=mu2, stringsAsFactors=FALSE);
		
		dftemp = dftemp[order(dftemp$time), ];
		dftemp$index = 1:nrow(dftemp);
		rownames(dftemp) = NULL;
		
		statevec = rep(1, nrow(phy$edge));

		if (nrow(dftemp) > 1)
		{
			for (k in 2:nrow(dftemp))
			{		
				s1 = which(phy$downseq == dftemp$node[k]);
				s2 = which(phy$downseq == phy$lastvisit[dftemp$node[k]]);
				descSet = phy$downseq[s1:s2];
				isDescendantNode = phy$edge[,2] %in% descSet;				
				statevec[isDescendantNode] = k;
			}				
		}

 		tmpEventSegMat = matrix(0, nrow=(max(phy$edge) + nrow(dftemp) - 2), ncol=4);
 		
		non.root = c(1:length(phy$tip.label), (length(phy$tip.label)+2):max(phy$edge));
		pos = 1;	
		
		is_noEventBranch = ! (phy$edge[,2] %in% dftemp$node);
		
		tmpEventSegMat[1:sum(is_noEventBranch), 1] = phy$edge[,2][is_noEventBranch];
		tmpEventSegMat[1:sum(is_noEventBranch),2] = phy$begin[is_noEventBranch];
 		tmpEventSegMat[1:sum(is_noEventBranch),3] = phy$end[is_noEventBranch];
 		tmpEventSegMat[1:sum(is_noEventBranch),4] = statevec[is_noEventBranch];		
 		
		eventnodeset = intersect(non.root, dftemp$node);
		pos = 1 + sum(is_noEventBranch);
		for (k in eventnodeset)
		{		
			events.on.branch = dftemp[dftemp$node == k, ];
			events.on.branch = events.on.branch[order(events.on.branch$time), ];
				
			fBranch = phy$edge[,2] == k;
 			start.times = c(phy$begin[fBranch], events.on.branch$time);
			stop.times = c(events.on.branch$time, phy$end[fBranch]);
			parent = phy$edge[,1][phy$edge[,2] == k];
			if (parent == (length(phy$tip.label) + 1))
			{
				# Parent is root:
				proc.set = c(1, events.on.branch$index);	
			}
			else
			{
				proc.set = c(statevec[phy$edge[,2] == parent], events.on.branch$index);			
			}
				
 			zzindex = pos:(pos+nrow(events.on.branch));	
				
			tmpEventSegMat[zzindex, 1] = rep(k, length(zzindex));
			tmpEventSegMat[zzindex, 2] = start.times;
			tmpEventSegMat[zzindex, 3] = stop.times;
			tmpEventSegMat[zzindex, 4] = proc.set;		
			pos = pos + 1 + nrow(events.on.branch);
		}
		
 		tmpEventSegMat = tmpEventSegMat[order(tmpEventSegMat[,1]), ];
 	
 		eventBranchSegs[[i]] = tmpEventSegMat;

		tipstates = numeric(length(phy$tip.label));
		tipstates = statevec[phy$edge[,2] <= phy$Nnode+1];
		tipstates = tipstates[order(phy$edge[phy$edge[,2] <= phy$Nnode+1,2])];
		
 		### Compute tip rates:
 
		stoptime = max(branching.times(phy));
		
		tiplam = dftemp$lam1[tipstates] * exp(dftemp$lam2[tipstates] * (stoptime - dftemp$time[tipstates]));
		tipmu = dftemp$mu1[tipstates];
		
		meanTipMu = meanTipMu + tipmu/nsamples;
		meanTipLambda = meanTipLambda + tiplam/nsamples;
		
		### List assignments and metadata across all events:
		eventData[[i]] = dftemp;	
		eventVectors[[i]]  = statevec;
		numberEvents[i] = nrow(dftemp);
		tipStates[[i]] = tipstates;
		
		tipLambda[[i]] = tiplam;
		tipMu[[i]] = tipmu;	
 	}
 	
	phy$numberEvents = numberEvents;
	phy$eventData = eventData;
	phy$eventVectors = eventVectors;
	phy$tipStates = tipStates;
	phy$tipLambda = tipLambda;
	phy$meanTipLambda = meanTipLambda;
	phy$eventBranchSegs = eventBranchSegs; 	
	phy$tipMu = tipMu;
	phy$meanTipMu = meanTipMu;
	if(type == 'diversification')
	{	
		phy$type = 'diversification';
	}
	else
	{
		phy$type = 'trait';	
	}
 	
	# Inherits attributes of class phylo
	# plus adds new class: bamm-data
	class(phy) = 'bammdata';
	return(phy);
}
