######################################
#	Internal function called by dtRates(...)
#
#
decimals = function(x)
{
	if(x%%1 != 0)
	{
		return(nchar(strsplit(as.character(x),".",fixed=TRUE)[[1]][[2]]));
	}	
	else
	{
		return(10);
	}
}