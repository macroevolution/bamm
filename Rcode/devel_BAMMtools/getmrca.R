getmrca = function(phy,t1,t2)
{
	ne = as.integer(dim(phy$edge)[1]);
	anc = as.integer(phy$edge[,1]);
	desc = as.integer(phy$edge[,2]);
	root = as.integer(length(phy$tip.label) + 1);
	
	.C('fetchmrca',anc,desc,root,ne,as.integer(t1),as.integer(t2),integer(1))[[7]];
}