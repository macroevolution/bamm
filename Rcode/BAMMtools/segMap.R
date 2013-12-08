######################################
#	Internal function called by dtRates(...)
#
#
segMap = function(nodes,begin,end,tau)
{
	foo = function(x,tau)
	{
		len = (x[3] - x[2])/tau; if (len%%1 == 0) len = len+1;
		ret = seq(x[2],x[3],length.out=len);
		if(length(ret) == 1) return(matrix(x,nrow=1));
		#ret = seq(x[2],x[3],length.out=length(ret));
		ret = rep(ret,each=2); ret=ret[2:(length(ret)-1)];
		ret = matrix(ret,ncol=2,byrow=TRUE);
		return(cbind(matrix(rep(as.integer(x[1]),nrow(ret)),ncol=1), ret));
	}
	times = cbind(nodes,begin,end);
	ret = apply(times,1,foo,tau);
	return(do.call(rbind,ret));	
}