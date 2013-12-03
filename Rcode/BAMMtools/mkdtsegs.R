######################################
#	Internal function called by plot.dtrates(...)
#	
#
mkdtsegs = function(x,tau)
{
	bn = sqrt((x[3]-x[1])^2 + (x[4]-x[2])^2);
	len = bn/tau; if (len %% 1 == 0) len = len + 1;
	
	j = seq(x[1],x[3],length.out=len);
	if(length(j) == 1) return(matrix(x[c(1,3,2,4,5)],nrow=1));
		
	k = seq(x[2],x[4],length.out = len);
	
	j = rep(j,each=2); j = j[2:(length(j)-1)];
	j = matrix(j,ncol=2,byrow=TRUE);	
	k = rep(k,each=2); k = k[2:(length(k)-1)];
	k = matrix(k,ncol=2,byrow=TRUE);	
	l = matrix(rep(x[5],nrow(j)),ncol=1);
	return(cbind(j,k,l));	
}