print.bammdata = function(ephy)
{
	print.phylo(as.phylo.bammdata(ephy));
	nsamples = length(ephy$eventData);
	cat(paste("\nPosterior samples:", nsamples,"\n\n"));
	cat("List elements:\n");
	cat("\t",names(ephy)[1:10]);
	cat("\t",names(ephy)[11:length(ephy)]);
}
