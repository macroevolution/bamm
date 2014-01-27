
library(BAMMtools);
data(whales, events.whales);
ewhales <- getEventData(whales, events.whales, burnin=0.1);
cstw <- cumulativeShiftProbsTree(ewhales);

ec <- rep('gray40', length(whales$edge.length));
ec[cstw$edge.length >= 0.95] <- 'red';

plot(whales, edge.color=ec, show.tip.label=F);




data(primates, events.primates);
eprimates <- getEventData(primates, events.primates, burnin=0.1, type = 'trait');
cstp <- cumulativeShiftProbsTree(eprimates);

ec <- rep('gray40', length(primates$edge.length));
ec[cstp$edge.length >= 0.95] <- 'red';

plot(primates, edge.color=ec, show.tip.label=F);

##################
data(prior.primates, mcmc.primates);
bf <- computeBayesFactors(mcmc.primates, prior.primates)

ec <- rep('gray40', length(primates$edge.length));
ec[cstp$edge.length >= 0.95] <- 'red';

plot(primates, edge.color=ec, show.tip.label=F);





