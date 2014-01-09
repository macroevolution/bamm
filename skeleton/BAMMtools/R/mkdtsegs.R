######################################
#	Internal function called by plot.dtrates(...)
#	
#
mkdtsegs = function(x,tau,phy,tH)
{
	#bn = sqrt((x[3]-x[1])^2 + (x[4]-x[2])^2);
	#len = bn/tau; if (len %% 1 == 0) len = len + 1;
	len = (phy$end[match(x[5],phy$edge[,2])]/tH-phy$begin[match(x[5],phy$edge[,2])]/tH)/tau; if (len %% 1 == 0) len = len + 1;
	
	j = seq(x[1],x[3],length.out=len);
	if(length(j) == 1) return(matrix(x[c(1,3,2,4,5)],nrow=1));
		
	k = seq(x[2],x[4],length.out = len);
	
	j = rep(j,each=2); j = j[-c(1,length(j))];
	j = matrix(j,ncol=2,byrow=TRUE);	
	k = rep(k,each=2); k = k[-c(1,length(k))];
	k = matrix(k,ncol=2,byrow=TRUE);	
	l = matrix(rep(x[5],nrow(j)),ncol=1);
	return(cbind(j,k,l));	
}