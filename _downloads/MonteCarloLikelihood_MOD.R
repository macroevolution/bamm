
MonteCarloLikelihood_MOD <- function(tree,
                                 sample,
                                 transition_rate,
                                 lambda_prior,
                                 mu_prior,
                                 num.intervals = 5000,
                                 reps = 50000) {
  
  mrcas <- mrca(tree)
  times <- branching.times(tree)
  
  # Get the processes
  num_processes <- nrow(sample)
  processes <- do.call(rbind,apply(sample,1,function(this_process){
    
    # Get the time of the process
    event_time <- max(times) - as.numeric(this_process[[4]])
    
    # Get the node immediately descending from the process
    if( event_time != max(times) ) {
      node <- mrcas[this_process[[2]],this_process[[3]]]
    } else {
      node <- tree$Nnode + 2
    }
    
    # Get the age of the node immediately descending from the process
    node_time <- as.numeric(times[as.character(node)])
    
    # Get the event parameters
    speciation <- as.numeric(this_process[[5]])
    extinction <- as.numeric(this_process[[7]])
    
    # Get the extinction probability
    simulated_times <- replicate(reps,SimulateCPBDP(time=event_time,init_lambda=speciation,init_mu=extinction,
                                                    transition_rate=transition_rate,lambda_prior=lambda_prior,
                                                    mu_prior=mu_prior))
    
    # Get the constant-rate extinction probability
    simulated_times_constant <- replicate(reps,SimulateCPBDP(time=event_time,init_lambda=speciation,init_mu=extinction,
                                                             transition_rate=0,lambda_prior=lambda_prior,
                                                             mu_prior=mu_prior))
    
    res <- data.frame(event_time=event_time,node=node,node_time=node_time,speciation=speciation,
                      extinction=extinction,
                      simulated_times=I(list(simulated_times)),
                      simulated_times_constant=I(list(simulated_times_constant)))
    
    return(res)
    
  }))
  
  # Now set up the likelihood computation
  interval_size      <- max(times) / (num.intervals - 1)
  intervals          <- seq(0,max(times),by=interval_size)
  post_order_nodes   <- as.numeric(names(sort(times,decreasing=FALSE)))
  d_nodes            <- rep(0,nrow(tree$edge))
  
  x <- length(post_order_nodes)
  y <- 0
  
  # Start all the tip probabilities at 1
  d_nodes[!tree$edge[,2] %in% tree$edge[,1]] <- 1
  
  for (this_node in post_order_nodes) {
    
    y <- y + 1
    
    which_edges        <- which(tree$edge[,1] %in% this_node)
    descendants        <- tree$edge[which_edges,]
    left_desc          <- descendants[1,2]
    right_desc         <- descendants[2,2]
    start_time         <- as.numeric(times[as.character(this_node)])
    
    if ( left_desc <= tree$Nnode + 1 ) {
      left_end_time <- 0.0
    } else {
      left_end_time <- as.numeric(times[as.character(left_desc)])
    }
    if ( right_desc <= tree$Nnode + 1 ) {
      right_end_time <- 0.0
    } else {
      right_end_time <- as.numeric(times[as.character(right_desc)])
    }
    
    left_intervals     <- getBranchIntervals(this_node,left_desc,start_time,left_end_time,num.intervals,processes,tree)
    left_d             <- d_nodes[which_edges[1]]
    left_d             <- computeBranchProbability(left_d,left_intervals,processes)
    
    right_intervals    <- getBranchIntervals(this_node,right_desc,start_time,right_end_time,num.intervals,processes,tree)
    right_d            <- d_nodes[which_edges[2]]
    right_d            <- computeBranchProbability(right_d,right_intervals,processes)
    
    new_d              <- left_d * right_d
    
    
    res <- list()
    
    if ( this_node != tree$Nnode + 2 ) {
      
      # Internal node
      node_process <- processes[getProcessForBranchAtTime(this_node,this_node,start_time,processes,tree),]
      d_nodes[tree$edge[,2] %in% this_node] <- new_d * node_process$speciation
      
    } else {
      
      # Root
      node_process <- processes[getProcessForBranchAtTime(this_node,this_node,start_time,processes,tree),]
      
      res$node_process <- node_process
      res$start_time <- start_time
      res$root_e_prob <-  getExtinctionProbability(node_process,start_time)
     
      
      log_likelihood <- log(new_d / (1 - getExtinctionProbability(node_process,start_time))^2)
      res$loglik <- log_likelihood
    }
    
  }
  
  return(res)
  
}
