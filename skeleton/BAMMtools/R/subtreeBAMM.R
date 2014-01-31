subtreeBAMM<-function(ephy,tips=NULL,node=NULL)
{
  if (! 'bammdata' %in% class(ephy)) {
    stop('Object phy must be of class bammdata');
  }
  if(is.null(tips)& is.null(node)){
    stop("need to specify either the tips or a innernode on the tree for subsetting");
  }
  if((!is.null(tips))&(! is.null(node))){
    cat("specified both tips and node, will subset the bammdata object according to the node\n")
  }
  # get the whole tree
  wtree <- as.phylo.bammdata(ephy)
  
  #wtree<-list()
  #wtree$edge<-phy$edge
  #wtree$Nnode<-phy$Nnode
  #wtree$tip.label<-phy$tip.label
  #wtree$edge.length<-phy$edge.length
  #class(wtree)<-'phylo'
  
  if(! is.null(node)){
    if(node<=length(wtree$tip.label)){
      stop('error: the node corresponds to one single tip')
    }else if (node> max(wtree$edge)){
      stop('the node does not exist on the tree')
    }
    
    tips <- ephy$downseq[which(ephy$downseq == node):which(ephy$downseq == ephy$lastvisit[node])]
    tips <- wtree$tip.label[tips[tips <= wtree$Nnode + 1]]
    #get_all_kid<-function(phy,node){
    #  if(node<=phy$Nnode+1){
    #    return(node);
    #  }else{
    #    tip<-c();
    #    innernode<-c(node);
    #    while (length(innernode)>0){
    #      l<-phy$edge[,1] %in% innernode
    #      m<-phy$edge[l,2]
    #      tip<-c(tip,m[m<=phy$Nnode+1])
    #      innernode<-m[m>phy$Nnode+1]
    #    }
    #    return(tip)
    #  }
    #}
    #tips<-wtree$tip.label[get_all_kid(wtree,node)]
  
  }else{
    if( length(tips)==1){
      stop("need more than one tip to subset the bammdata object")
    }
    if(class(tips) != 'character'){
      #stop("tips should be a character string")
      tips <- wtree$tip.label[tips]
    }
  }
  
  #get the sub tree
  stree<-drop.tip(wtree,tip=wtree$tip.label[! wtree$tip.label %in% tips])
  sNtip<-length(stree$tip.label)
  #for every node on the subset tree, get its corresponding node on the whole tree
  oldnode<-integer(length=max(stree$edge))
  oldnode[1:sNtip]<-match(stree$tip.label,wtree$tip.label)
  # for every inner node on the stree get its left kid and right kid
  lrkids<-matrix(0,nrow=0,ncol=3);  
  get_lr_kids<-function(tree,node){
    ll<-which(tree$edge[,1]==node)
    if(length(ll)>0){
      l<-get_lr_kids(tree,tree$edge[ll[1],2])
      r<-get_lr_kids(tree,tree$edge[ll[2],2])
      return(rbind(l,r,c(node,min(l),min(r))))
    }else{
      return(matrix(c(node,node,node),nrow=1,ncol=3))
    }
  }
  lrkids<-get_lr_kids(stree,sNtip+1);
  lrkids<-lrkids[which(lrkids[,1]>sNtip),]
  oldnode[lrkids[,1]]<-getmrca(wtree,oldnode[lrkids[,2]],oldnode[lrkids[,3]])
  # for every branch on the new tree, get its corresponding branches on the whole tree
  # innernodes on the original trees shift towards the tip
  inner_node_shift<-matrix(0L,nrow=0,2)
  for (x in 1: dim(stree$edge)[1]){
    downnode<-oldnode[stree$edge[x,2]]
    upnode<-oldnode[stree$edge[x,1]]
    j<-downnode;
    j=wtree$edge[which(wtree$edge[,2]==j),1]
    while (j !=upnode){
      inner_node_shift<-rbind(inner_node_shift,c(j,stree$edge[x,2]))
      j=wtree$edge[which(wtree$edge[,2]==j),1]
      #cat(j)
    }
  }
  
  subphy<-stree;
  class(subphy)<-"bammdata"
  stree<-getStartStopTimes(stree)
  subphy$begin<-stree$begin
  subphy$end<-stree$end
  stree <- getRecursiveSequence(stree);
  subphy$downseq<-stree$downseq
  subphy$lastvisit<-stree$lastvisit
  subphy$numberEvents<-c()
  subphy$eventData<-list()
  subphy$eventVectors<-list()
  subphy$tipStates<-list()
  subphy$eventBranchSegs<-list()
  subphy$tipLambda<-sapply(ephy$tipLambda,function(x){x[oldnode[1:sNtip]]})
  subphy$meanTipLambda<-ephy$meanTipLambda[oldnode[1:sNtip]]
  subphy$tipMu<-sapply(ephy$tipMu,function(x){x[oldnode[1:sNtip]]})
  subphy$meanTipMu<-ephy$meanTipMu[oldnode[1:sNtip]]
  subphy$type<-ephy$type
  
  
  if(oldnode[sNtip+1]==length(ephy$tip.label)+1){
    rootshift<-0
  }else{
    rootshift<-max(ephy$end)-max(subphy$end)
  }
  for (en in 1:length(ephy$numberEvents)){
    #cat(en,"\n")
    eventData<-ephy$eventData[[en]];
    eventVectors<-ephy$eventVectors[[en]]
    if(rootshift>0){
      root_process<-eventVectors[which(ephy$edge[,2]==oldnode[sNtip+1])]
    }else{
      root_process<-1L
    }
    eventData$node<-sapply(eventData$node,function(x){
      if(x %in% oldnode){
        which(oldnode==x);
      }else if(x %in% inner_node_shift[,1]){
        inner_node_shift[ which(inner_node_shift[,1]==x),2]
      }else{
        0L
      }
    })
    eventData$time<-eventData$time-rootshift
    eventData<-eventData[(eventData$node>0)|(eventData$index==root_process),]
    #fix the root process
    eventData$node[eventData$node==0L]<-sNtip+1; #the only reason that a event that does not happen on stree and still is kept is that it is a root process
    l<-which(eventData$time<0);
    if(length(l)>0){
      if(length(l)>1){
        l<-which(eventData$time==max(eventData$time[l]))
      }
      #eventData$lam1[l]<-eventData$lam1[l]*exp((0-eventData$time[l])*eventData$lam2[l])
      eventData$lam1[l]<-exponentialRate(0-eventData$time[l], eventData$lam1[l], eventData$lam2[l])
      eventData$time[l]<-0;
    }
    eventData<-eventData[eventData$time>=0,]
    eventData<-eventData[order(eventData$time),]
    
    #now fix the index of the process
    subphy$numberEvents[en]<-dim(eventData)[1]
    newprocess<-eventData$index;
    eventData$index<-1:dim(eventData)[1];
    subphy$eventData[[en]]<-eventData
    #fix the tipstate
    tipState<-ephy$tipStates[[en]];
    tipState<-tipState[oldnode[1:sNtip]]
    subphy$tipState[[en]]<-match(tipState,newprocess)
    #fix the eventVectors
    eventVectors<-sapply(subphy$edge[,2], function(x){eventVectors[which(ephy$edge[,2]==oldnode[x])]})
    subphy$eventVectors[[en]]<-match(eventVectors,newprocess)
    #write up a new eventBranchSeg
    eventBranchSegs<-cbind(subphy$edge[,2],subphy$begin,subphy$end,subphy$eventVectors[[en]])
    eventBranchSegs<-eventBranchSegs[ ! eventBranchSegs[,1] %in% eventData$node[eventData$time>0],]
    if(length(eventData$node[eventData$time>0])>0){
      for (n in unique(eventData$node[eventData$time>0])){
        branch_event<-eventData[eventData$node==n,]
        branch_event<-branch_event[order(branch_event$time),]
        upnode<-subphy$edge[ which(subphy$edge[,2]==n),1]
        if(upnode==sNtip+1){startproc<-1L;}else{startproc<-subphy$eventVectors[[en]][which(subphy$edge[,2]==upnode)]}
        starttime<-subphy$begin[which(subphy$edge[,2]==n)];
        for (be in 1:dim(branch_event)[1]){
          #l<-c(n,starttime,branch_event$time[be],startproc);
          eventBranchSegs<-rbind(eventBranchSegs,c(n,starttime,branch_event$time[be],startproc));
          startproc<-branch_event$index[be]
          starttime<-branch_event$time[be]
        }    
        eventBranchSegs<-rbind(eventBranchSegs,c(n,starttime,subphy$end[which(subphy$edge[,2]==n)],startproc));
      }
    }
    subphy$eventBranchSegs[[en]]<-eventBranchSegs[order(eventBranchSegs[,1]),]
  }
  return(subphy);
}
