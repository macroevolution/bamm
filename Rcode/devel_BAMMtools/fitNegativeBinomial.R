
# This is a quick fix. Needs to be made more robust to optimization issues. 
# In particular, this will probably not work on large datasets 
#	eg if someone feeds in a vector
# 	of 1000000 observations
#
#  This is a PRIVATE function in BAMMtools, used 
#	only by computeBayesFactors

fitNegativeBinomial <- function(x){

	fx <- function(p){
		sz <- exp(p[1]);
		mu <- exp(p[2]);
		return(-sum(dnbinom(x, size=sz, mu=mu, log=T)));
	}
	
	ss <- c(1, 1);
	ll <- 0;
	class(ll) <- 'try-error';
	while (class(ll) == 'try-error'){
		ss <- c(runif(1, 0, 100), runif(1, 0, 100));
		ll <- try(fx(log(ss)));
	}
	
	res <- optim(par=log(ss), fn=fx, method = 'Nelder');
	
	obj <- list();
	obj$logLik <- -res$value;
	obj$size <- exp(res$par[1]);
	obj$mu <- exp(res$par[2]);
	obj$conv <- res$convergence;
	return(obj);
	
}



