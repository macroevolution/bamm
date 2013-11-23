
wd <- getwd();
wd_split <- unlist(strsplit(wd, split='/'));

if (wd_split[length(wd_split)] != 'Rcode'){
	stop("This function only works from the bamm/Rcode directory\n");
} 

ll <- dir('BAMMtools/');

ll <- paste('BAMMtools/', ll, sep='');
for (i in ll){
	source(i);
}
 
 
rm(ll, i, wd, wd_split);