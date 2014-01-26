library(BAMMtools);

tmptree <- '((A:2,B:2):4,((C:1.5,D:1.5):1.5,(E:2,F:2):1):3);'

v <- read.tree(text=tmptree);

d1 <- getDesc(v, node=8)$desc_set;
d2 <- getDesc(v, node=9)$desc_set;

edgecolor <- rep('black', nrow(v$edge));

edgecolor[v$edge[,1] %in% d1] <- 'red';
edgecolor[v$edge[,1] %in% d2] <- 'blue';

quartz.options(height=6, width=6);
par(oma=c(0,0,0,0));
plot.phylo(v, edge.width=3, edge.color=edgecolor, cex=2, label.offset=0.15);











