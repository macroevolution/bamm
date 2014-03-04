speciesByRatesMatrix = function(ephy, nslices, index = NULL, spex = "s") {
	phy = as.phylo.bammdata(ephy);
	seq.nod = .Call("seq_root2tip", phy$edge, length(phy$tip.label), phy$Nnode, PACKAGE = "BAMMtools");
	if (nslices <= 100) {
		tvec = (seq(0, 1, 0.01)+0.005) * max(branching.times(phy));
		tvec = tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy = dtRates(ephy, 0.01, index, tmat = TRUE);
	}
	else if (nslices > 100 && nslices <= 500) {
		tvec = (seq(0, 1, 0.002)+0.001) * max(branching.times(phy));
		tvec = tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy = dtRates(ephy, 0.002, index, tmat = TRUE);
	}
	else if (nslices > 500 && nslices <= 1000) {
		tvec = (seq(0, 1, 0.001)+0.0005) * max(branching.times(phy));
		tvec = tvec[seq.int(1,length(tvec),length.out=nslices+1)];
		ephy = dtRates(ephy, 0.001, index, tmat = TRUE);
	}
	else {
		stop("Max slices (1000) exceeded.  Choose a smaller number of slices");
	}
	ret = lapply(seq.nod, function(x) {
		path = which(ephy$dtrates$tmat[,1] %in% x);
		ids = sapply(tvec[-length(tvec)], function(y) which(ephy$dtrates$tmat[path,2] <= y & ephy$dtrates$tmat[path,3] > y));
		if (ephy$type == "trait") {
			ephy$dtrates$rates[path][ids];
		}
		else {
			if (tolower(spex) == "s") {
				ephy$dtrates$rates[[1]][path][ids];
			}
			else if (tolower(spex) == "e") {
				ephy$dtrates$rates[[2]][path][ids];
			}
			else {
				ephy$dtrates$rates[[1]][path][ids] - ephy$dtrates$rates[[2]][path][ids];
			}
		}
	});
	ret = do.call(rbind, ret);
	rownames(ret) = phy$tip.label;
	return(list(times = tvec[-length(tvec)],rates = ret));	
}