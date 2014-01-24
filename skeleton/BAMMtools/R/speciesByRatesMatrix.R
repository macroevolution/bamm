speciesByRatesMatrix = function(ephy, nslices, index = NULL, spex = "s") {
	phy = as.phylo.bammdata(ephy);
	seq.nod = .Call("seq_root2tip", phy$edge, length(phy$tip.label), phy$Nnode, PACKAGE = "ape");
	if (nslices > length(phy$tip.label)) {
		warning("You are generating more time slice variables than species. Consider choosing a smaller number of slices");
	}
	if (nslices <= 100) {
		tvec = (seq(0, 1, 0.01)+0.005) * max(branching.times(phy));
		tvec = tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		m = dtRates(ephy, 0.01, index, tmat = TRUE);
	}
	else if (nslices > 100 && nslices <= 500) {
		tvec = (seq(0, 1, 0.002)+0.001) * max(branching.times(phy));
		tvec = tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		m = dtRates(ephy, 0.002, index, tmat = TRUE);
	}
	else if (nslices > 500 && nslices <= 1000) {
		tvec = (seq(0, 1, 0.001)+0.0005) * max(branching.times(phy));
		tvec = tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		m = dtRates(ephy, 0.001, index, tmat = TRUE);
	}
	else {
		stop("Max slices (1000) exceeded.  Choose a smaller number of slices");
	}
	tol = 0.0001*max(branching.times(phy));
	ret = lapply(seq.nod, function(x) {
		path = which(m$dtrates$tmat[,1] %in% x);
		ids = unlist(sapply(tvec[-length(tvec)], function(y) which(m$dtrates$tmat[path,2] <= y & m$dtrates$tmat[path,3] > y)));
		if (ephy$type == "trait") {
			return(m$dtrates$rates[ids]);
		}
		else {
			if (tolower(spex) == "s") {
				return(m$dtrates$rates[[1]][ids]);
			}
			else if (tolower(spex) == "e") {
				return(m$dtrates$rates[[2]][ids]);
			}
			else {
				return(m$dtrates$rates[[1]][ids] - m$dtrates$rates[[2]][ids]);
			}
		}
	});
	ret = do.call(rbind, ret);
	rownames(ret) = phy$tip.label;
	return(ret);	
}