
require(geiger)
require(diversitree)


####-----------------------------------------
#### The following block of functions:
#
# simulateCBDPTtree
# CPBDPStochasticMap
# pruneCPBDPTree
# computeShiftDescendants
#
# were taken from the Dryad supplement to accompany Moore et al
# PNAS 2016: www.pnas.org/cgi/doi/10.1073/pnas.1518659113
# 
# Dryad doi: http://datadryad.org/resource/doi:10.5061/dryad.mb0sd
#
# We slightly modified SimulateCBDPTree such that it would 
#   generate a larger range of tree sizes by removing the 
#   size selection bias that contributed to ascertainment bias 
#   in the MEA PNAS article
 
SimulateCBDPTree <- function(time,
                             transition_rate,
                             lambda_function,
                             mu_function,
                             init_lambda = NULL,
                             init_mu = NULL,
                             verbose = TRUE,
                             condition_on_survival = TRUE,
                             condition_on_root = TRUE,
                             NMAX,
                             MAX_FAILS = 1,
                             cur_fails = 0) {

  # If initial lambda and mu are not provided, simulate them from the prior distribution
  if ( is.null(init_lambda) ) {
    initial_speciation <- lambda_function()
  } else {
    initial_speciation <- init_lambda
  }

  if ( is.null(init_mu) ) {
    initial_extinction <- mu_function()
  } else {
    initial_extinction <- init_mu
  }

  # Initialize the number of species
  lineage_index <- 3
  current_num_species <- 2
  total_number_of_transitions <- 0

  # Initialize the edge matrix
  edge <- data.frame(ancestor=1,descendant=2:3,start_time=0,end_time=NA,
                     current_speciation_rate=initial_speciation,current_extinction_rate=initial_extinction,
                     speciation_rates=I(list(initial_speciation,initial_speciation)),extinction_rates=I(list(initial_extinction,initial_extinction)),
                     transition_times=I(list(NA,NA)),status="alive",states=I(list(0,0)),stringsAsFactors=FALSE)

  # Start the simulation
  current_time <- 0

  if(verbose) bar <- txtProgressBar(style=3,width=40)
 

  fail_count = cur_fails

  while( TRUE ) {

    # get the next event time
    these_rates <- (edge$current_speciation_rate + edge$current_extinction_rate + transition_rate) * as.numeric(edge$status == "alive")
    next_times <- suppressWarnings(rexp(length(these_rates),these_rates))
    which_edge_has_event <- which.min(next_times)
    next_event_time <- next_times[which_edge_has_event]

    # increment the current time; if it is greater than max time, break
    current_time <- current_time + next_event_time

    if ( current_time > time ) {
      current_time <- time
      if(verbose) setTxtProgressBar(bar,current_time/time)
      break
    }

    if(verbose) setTxtProgressBar(bar,current_time/time)

    # perform the next event
    this_edge <- edge[which_edge_has_event,]
    these_rates <- c(this_edge$current_speciation_rate,this_edge$current_extinction_rate,transition_rate)
    total_rate <- sum(these_rates)
    probs <- these_rates / total_rate
    probs[these_rates == Inf] <- 1.0
    next_event_type <- sample.int(3,1,prob=probs)

    if ( next_event_type == 1 ) {

      # Speciation event
      edge[which_edge_has_event,]$status <- "speciated"
      edge[which_edge_has_event,]$end_time <- current_time

      this_ancenstor <- edge[which_edge_has_event,]$descendant
      these_states <- edge[which_edge_has_event,]$states[[1]]
      this_state <- these_states[length(these_states)]
      new_descendants <- 1:2 + lineage_index

      these_species <- data.frame(ancestor=this_ancenstor,descendant=new_descendants,start_time=current_time,end_time=NA,
                                  current_speciation_rate=these_rates[1],current_extinction_rate=these_rates[2],
                                  speciation_rates=I(list(these_rates[1],these_rates[1])),extinction_rates=I(list(these_rates[2],these_rates[2])),
                                  transition_times=I(list(NA,NA)),status="alive",states=I(list(this_state,this_state)),stringsAsFactors=FALSE)

      edge <- rbind(edge,these_species)

      current_num_species <- current_num_species + 1
      lineage_index <- lineage_index + 2

      # cat("current time: ",current_time," -- speciation event -- current species: ",current_num_species," -- number of transitions: ",total_number_of_transitions,"\n",sep="")

    } else if ( next_event_type == 2 ) {

      # Extinction event
      edge[which_edge_has_event,]$status <- "extinct"
      edge[which_edge_has_event,]$end_time <- current_time
      current_num_species <- current_num_species - 1

      # cat("current time: ",current_time," -- extinction event -- current species: ",current_num_species," -- number of transitions: ",total_number_of_transitions,"\n",sep="")

    } else if ( next_event_type == 3 ) {

      # Transition event
      total_number_of_transitions <- total_number_of_transitions + 1
      new_speciation_rate <- lambda_function()
      new_extinction_rate <- mu_function()

      these_speciation_rates <- this_edge$speciation_rates
      these_extinction_rates <- this_edge$extinction_rates
      these_transition_times <- this_edge$transition_times
      these_states           <- this_edge$states

      these_speciation_rates[[1]] <- c(these_speciation_rates[[1]],new_speciation_rate)
      these_extinction_rates[[1]] <- c(these_extinction_rates[[1]],new_extinction_rate)
      these_transition_times[[1]] <- c(these_transition_times[[1]],current_time)
      these_states[[1]]           <- c(these_states[[1]],total_number_of_transitions)

      this_edge$current_speciation_rate <- new_speciation_rate
      this_edge$current_extinction_rate <- new_extinction_rate
      this_edge$speciation_rates <- these_speciation_rates
      this_edge$extinction_rates <- these_extinction_rates
      this_edge$transition_times <- these_transition_times
      this_edge$states           <- these_states
      edge[which_edge_has_event,] <- this_edge

      # cat("current time: ",current_time," -- transition event -- current species: ",current_num_species," -- number of transitions: ",total_number_of_transitions,"\n",sep="")

    }

    # if the total number of species is 0, break
    if (current_num_species == 0) {
      if(verbose) {
        flush.console()
        if (verbose) cat("\nTree died!\n")
      }
   
      
      break
    }

    if ( current_num_species > NMAX ) {
      cat("Aborted / N = ", current_num_species, " Fails: ", fail_count, "\t" )
      
      fail_count <- fail_count + 1
      if (fail_count >= MAX_FAILS){
      	return(NA)
      }
      
      break
    }
 
  }

  if ( current_num_species < 1 & condition_on_survival ) {
    return(SimulateCBDPTree(time=time,transition_rate=transition_rate,
                             lambda_function=lambda_function,
                             mu_function=mu_function,
                             init_lambda=init_lambda,
                             init_mu=init_mu,
                             verbose=verbose,
                             condition_on_survival=condition_on_survival, NMAX=NMAX, MAX_FAILS = MAX_FAILS, cur_fails = fail_count))
  }

  if ( current_num_species > NMAX ) {
    return(SimulateCBDPTree(time=time,transition_rate=transition_rate,
                            lambda_function=lambda_function,
                            mu_function=mu_function,
                            init_lambda=init_lambda,
                            init_mu=init_mu,
                            verbose=verbose,
                            condition_on_survival=condition_on_survival, NMAX=NMAX, MAX_FAILS = MAX_FAILS, cur_fails = fail_count))
  }

  edge$end_time[is.na(edge$end_time)] <- current_time

  # Now make the tree ape-friendly
  phy_edge_length <- edge$end_time - edge$start_time
  phy_edge <- matrix(NA,nrow=nrow(edge),ncol=2)

  new_tips <- 1:sum(edge$status != "speciated")
  new_root <- max(new_tips) + 1
  new_nodes <- new_root + 1:(length(new_tips)-2)
  old_nodes <- edge$descendant[edge$status == "speciated"]

  phy_edge[,1] <- new_nodes[match(edge$ancestor,old_nodes)]
  phy_edge[,2] <- new_nodes[match(edge$descendant,old_nodes)]

  phy_edge[edge$status != "speciated",2] <- new_tips
  phy_edge[edge$ancestor == 1,1] <- new_root

  Nnode <- new_root - 2

  phy <- list()
  class(phy) <- "phylo"
  phy$edge <- phy_edge
  phy$edge.length <- phy_edge_length
  phy$tip.label <- new_tips
  phy$Nnode <- Nnode
  phy$full_process <- edge

  if ( condition_on_root ) {

    these_descendants <- phy$edge[phy$edge[,1] == Nnode + 2,2]
    left_descendants <- tips(phy,these_descendants[1])
    right_descendants <- tips(phy,these_descendants[2])
    if ( !any("alive" %in% edge[phy$edge[,2] %in% left_descendants,]$status) | !any("alive" %in% edge[phy$edge[,2] %in% right_descendants,]$status) ){
      if (verbose) cat("\nRoot died!\n")
      return(SimulateCBDPTree(time=time,transition_rate=transition_rate,
                              lambda_function=lambda_function,
                              mu_function=mu_function,
                              init_lambda=init_lambda,
                              init_mu=init_mu,
                              verbose=verbose,
                              condition_on_survival=condition_on_survival, NMAX = NMAX, MAX_FAILS = MAX_FAILS, cur_fails = fail_count))
    }
  }

  return(phy)

}

CPBDPStochasticMap <- function(tree) {

  edge <- tree$full_process

  # Now make the SIMMAP-formated transition events
  total_states <- 1 + max(unlist(edge$states))

  maps <- lapply(1:nrow(edge),function(x){
    row <- edge[x,]
    transition_times <- row$transition_times[[1]]
    durations <- diff(na.omit(c(row$start_time,transition_times,row$end_time)))
    names(durations) <- row$states[[1]]
    return(durations)
  })

  blank_map <- rep(0,total_states)
  names(blank_map) <- 1:total_states-1

  mapped.edge <- do.call(rbind,lapply(maps,function(map){
    for(i in 1:length(map)) blank_map[names(map[i])] <- blank_map[names(map[i])] + map[i]
    return(blank_map)
  }))

  tree$maps <- maps
  tree$mapped.edge <- mapped.edge

  return(tree)

}

pruneCPBDPTree <- function(tree) {

  edge <- tree$full_process

  # these nodes are already sorted from youngest to oldest
  # we reverse them to get the downpass sequence
  # exclude the root
  nodes <- rev(unique(edge$ancestor)[-1])

  for ( node in nodes ) {

    node_ancestor <- edge$ancestor[edge$descendant == node]
    these_edges <- edge[edge$ancestor %in% node,]
    if ( all(these_edges$status == "extinct") ) {
      # if both of the nodes descendants are extinct, drop all of them
      # and set the current node's status to extinct
      edge[edge$descendant %in% node,]$status <- "extinct"
      edge <- edge[!edge$ancestor %in% node,]
    } else if ( any(these_edges$status == "extinct") ) {
      # otherwise, just drop the extinct tip
      extinct_edge <- these_edges[these_edges$status == "extinct",]
      extinct_descendant <- extinct_edge$descendant

      # remove the extinct descendant
      edge <- edge[edge$descendant != extinct_descendant,]

      # relabel the ancestor and descendants
      alive_edge <- these_edges[these_edges$status != "extinct",]
      alive_descendant <- alive_edge$descendant
      which_recipient <- which(edge$descendant == node)
      which_donor <- which(edge$descendant == alive_descendant)
      edge[which_recipient,]$end_time <- edge[which_donor,]$end_time
      edge[which_recipient,]$status <- edge[which_donor,]$status

      if ( length(na.omit(edge[which_donor,]$transition_times[[1]])) > 0 ) {
        # if the donor branch has transitions, glue together the rates
        edge[which_recipient,]$speciation_rates <- list(c(edge[which_recipient,]$speciation_rates[[1]],edge[which_donor,]$speciation_rates[[1]][!is.na(edge[which_donor,]$transition_times[[1]])]))
        edge[which_recipient,]$extinction_rates <- list(c(edge[which_recipient,]$extinction_rates[[1]],edge[which_donor,]$extinction_rates[[1]][!is.na(edge[which_donor,]$transition_times[[1]])]))
        edge[which_recipient,]$current_speciation_rate <- edge[which_donor,]$current_speciation_rate
        edge[which_recipient,]$current_extinction_rate <- edge[which_donor,]$current_extinction_rate
        edge[which_recipient,]$states <- list(c(edge[which_recipient,]$states[[1]],edge[which_donor,]$states[[1]][-1]))
        edge[which_recipient,]$transition_times <- list(as.numeric(c(edge[which_recipient,]$transition_times[[1]],na.omit(edge[which_donor,]$transition_times[[1]]))))
      }

      # now, remove the donor
      edge <- edge[-which_donor,]
      edge$descendant[edge$descendant == node] <- alive_descendant

    }

  }

  # remember to do something with the root
  node <- 1
  these_edges <- edge[edge$ancestor %in% node,]

  if ( all(these_edges$status == "extinct") ) {
    message <- "The whole tree went extint. What happened??"
    cat(message)
    class(message) <- "try-error"
    return(invisible(message))
  } else if ( any(these_edges$status == "extinct" ) ) {
    # We need to remove the root
    warning("Removing root. Was this supposed to happen?")
    edge <- edge[!edge$ancestor %in% node,]
  }

  # Now make the tree ape-friendly
  smallest_node <- min(edge$ancestor)
  phy_edge_length <- edge$end_time - edge$start_time
  phy_edge <- matrix(NA,nrow=nrow(edge),ncol=2)
  new_tips <- 1:sum(edge$status != "speciated")
  new_root <- max(new_tips) + 1
  new_nodes <- new_root + 1:(length(new_tips)-2)
  old_nodes <- edge$descendant[edge$status == "speciated"]
  phy_edge[,1] <- new_nodes[match(edge$ancestor,old_nodes)]
  phy_edge[,2] <- new_nodes[match(edge$descendant,old_nodes)]
  phy_edge[edge$status != "speciated",2] <- new_tips
  phy_edge[edge$ancestor == smallest_node,1] <- new_root
  Nnode <- new_root - 2

  phy <- list()
  class(phy) <- "phylo"
  phy$edge <- phy_edge
  phy$tip.label <- new_tips
  phy$edge.length <- phy_edge_length
  phy$Nnode <- Nnode
  phy$full_process <- edge

  if(!is.null(tree$maps)){
    phy <- CPBDPStochasticMap(phy)
  }

  return(phy)


}

computeShiftDescendants <- function(tree) {

  edge <- tree$full_process
  edge$extant <- apply(edge,1,function(row){
    descendants <- tips(tree,tree$edge[edge$descendant == row$descendant,2])
    any(edge[tree$edge[,2] %in% descendants,]$status == "alive")
  })

  tree$full_process <- edge

  return(tree)

}


#### End MEA block of functions
####-----------------------------------------


#### End block of BAMM-testing code
####-----------------------------------------


# Function to extract basic shift regime data
# including numbers of tips in each rate regime
# and number of branches included in each rate regime
# from a bammdata object.
#  This only extracts from a single indexed shift 
#  configuration. 
shiftRegimeData <- function(ed, index = 1){
	
	ex <- ed$eventData[[index]]
	tipstates <- table(ed$tipStates[[index]])
	tipstates <- tipstates[as.character(ex$index)]
	
	edgestates <- table(ed$eventVectors[[1]])
	edgestates <- edgestates[as.character(ex$index)]
	
	ex$ntips      <- tipstates
	ex$nbranches  <- edgestates 
	return(ex)
}


colorByRegime <- function(ed){
	
	ev <- ed$eventVectors[[1]]
	cols <- sample(rainbow(n = length(unique(ev)) ))
	
	colset <- cols[ev]
	
	return(colset)
}


birthdeath_fit <- function(phy, node = NULL, dropnodes = NULL){
	
	droptips <- NULL
	
	if (!is.null(dropnodes)){
		
		for (i in 1:length(dropnodes)){
			
			if (i <= length(phy$tip.label)){
				droptips <- c(droptips, phy$tip.label[i])
			}else{
				tmp <- extract.clade(phy, node = dropnodes[i])$tip.label
				droptips <- c(droptips, tmp)
			}
			
			
		}
		
		phy <- drop.tip(phy, droptips)
		
	}
	if (!is.null(node)){
	     phy <- extract.clade(phy, node = node)	
	}
	
	lfx <- make.bd(phy)
	lfunc <- function(par){
		 -lfx(exp(par))
	} 
	
	
	ntips <- length(phy$tip.label)
	
	init <- 1.25 * (log(ntips) - log(2)) / max(branching.times(phy))
	init <- log(c(init, 0.5*init))
	names(init) <- c("lambda", "mu")
	res <- optim(par=init, lfunc, method="Nelder" )
	
	rates <- round(exp(res$par), 3)
	names(rates) <- c("lambda", "mu")
	
	return(rates)
	
}



plotSetup <- function(mx = 0.8){
	plot.new()
	par(oma=c(1,1,1,1))
	par(mar=c(5,5,1,1))
	plot.window(xlim=c(0, mx), ylim=c(0, mx), asp=1)	
	abline(0, 1, lwd=3, col="gray50")
	axis(1, at=seq(-0.2, 1.2, by=0.2))
	axis(2, at=seq(-0.2, 1.2, by=0.2), las=1)
	mtext("True rate", side=1, line=3, cex=1.4)
	mtext("BAMM estimated rate", side = 2, line = 3, cex = 1.4)
	
	lines(x=c(0,0.1), y=c(0.6, 0.6), lwd=3, col="gray50")
	text(x=0.1, y=0.6, pos=4, label = "1:1 line", font = 3)
	lines(x=c(0,0.1), y=c(0.55, 0.55), lty="dotted", lwd=3)
	text(x=0.1, y=0.55, pos=4, label = "BAMM fit", font=3)
		
	
}

# There is a huge amount of 
#  overplotting in the following plots, 
#   since most taxa have identical rates.
#  So we will "jitter" the points slightly so you can see
#  roughly how many branches fall out together

j <- function(x){
	return(jitter(x, amount = 0.015))
}
 


# This function converts the events from a tree simulated using the MEA R simulator to BAMM's event data format
mea_to_edata <- function(tree)	{
	require(phytools)
	N <- Ntip(tree)
	SpecRates <- tree$full_process$current_speciation_rate
	tipLambda <- setNames(SpecRates[sapply(1:N, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)

	ExtRates <- tree$full_process$current_extinction_rate
	tipMu <- setNames(ExtRates[sapply(1:N, function(x,y) which(y==x), y=tree$edge[,2])], tree$tip.label)
	
	Times <- sapply(tree$full_process$transition_times, function(x) x[2])
	
	eventData <- data.frame(generation=0, leftchild="drop", rightchild="drop", abstime=-1, lambdainit=0.1, lambdashift=0, muinit=0.1, mushift=0, n_taxa=0)
	for (regime in unique(tipLambda))	{
		Include <- which(tipLambda == regime)
		Tips <- names(tipLambda)[Include]
		abstime <- Times[which(SpecRates==regime)[1]]
		if (is.na(abstime))	{
			abstime <- 0
		}
		if (length(Tips) > 1)	{
			Node <- findMRCA(tree, Tips)
			Descs <- getDescendants(tree, Node)
			Descs <- Descs[Descs <= N]
			lchild <- as.character(Descs[1])
			rchild <- as.character(Descs[length(Descs)])
			eventData <- rbind(eventData, data.frame(generation=0, leftchild=lchild, rightchild=rchild, abstime=abstime, lambdainit=regime, lambdashift=0, muinit=tipMu[lchild], mushift=0, n_taxa=length(Tips)))
		}
		else	 if (length(Tips) == 1){
			eventData <- rbind(eventData, data.frame(generation=0, leftchild=as.character(Tips), rightchild="NA", abstime=abstime, lambdainit=regime, lambdashift=0, muinit=tipMu[Tips], mushift=0, n_taxa=1))
		}		
	}

	eventData <- eventData[-1,]
	eventData <- eventData[order(eventData$abstime),]
	return(eventData)
}











