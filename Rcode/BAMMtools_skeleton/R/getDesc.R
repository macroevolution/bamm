getDesc <-
function(phy, node)
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
