getmrca = function(phy,t1,t2)
{
	ne = as.integer(dim(phy$edge)[1]);
	anc = as.integer(phy$edge[,1]);
	desc = as.integer(phy$edge[,2]);
	
	.C('fetchmrca',anc,desc,ne,as.integer(t1),as.integer(t2),integer(1))[[6]];
}