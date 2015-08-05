tn_state <- function(data,time,tr,unit){
  window <- c()
  begin <- F
  end <- F
  for(i in 1:length(data[,3])){
    t <- as.POSIXlt(data[i,3])
    diff <- as.numeric(difftime(time,t,units=unit))
    if(t<time && diff<tr){
      begin <- T
      window <- c(window,i)
    }else{
      if(begin) break
    }
  }
  window
}


findAllNodes <- function(data){
  tn <- graph.empty(directed=FALSE)
  for(i in 1:nrow(data)){
    s = toString(data[i,1])
    r = toString(data[i,2])
    #add new nodes if not exist
    if(length(which(V(tn)$name == s))==0){
      tn <- tn + vertex(s, infected=FALSE)
      print(paste("Node",s,"added.",sep=" "))
    }    
    if(length(which(V(tn)$name == r))==0){
      tn <- tn + vertex(r, infected=FALSE)
      print(paste("Node",r,"added.",sep=" "))
    }
  }
  return(tn)
}


recomputeWeights <- function(n, dt, norm, tr){
  todelete <- list()
  for(edge in E(n)){
    nlt <- E(n)[edge]$lt + (dt/norm)
    nw <- exp(-nlt)
    E(n)[edge]$weight <- nw
    if(nw < tr){
      todelete <- append(todelete,edge)
    }else{
      E(n)[edge]$lt <- nlt
    }
  }  
  n <- delete.edges(n,todelete)
  return(n)
}


countWeights <- function(n,bins){
  if(ecount(n)>0){
    h=hist(E(n)$weight,breaks=bins,plot=FALSE)
    return(h$counts)
  }else{
    return(matrix(0,ncol=bins))
  }
}


countActiveNodes <- function (n){
  an = 0
  components = decompose.graph(n)
  cn = 0
  for(x in components){
    if(ecount(x)>0){
      an = an + vcount(x)
      cn = cn + 1
    }
  }
  return(c(an,cn))
}


runEvents <- function(tr, norm, unit, bins, window, data){
  lastEventTime = as.POSIXct(data[1,3])
  anCnPerEvent = matrix(c(0,0),ncol=2)
  weightHist = matrix(ncol = length(bins))
  trans = c()
  for(i in 1:nrow(data[window,])){ 
    
    #read next event    
    s = toString(data[i,1])
    r = toString(data[i,2])
    t = as.POSIXct(data[i,3])
    
    #time since last event
    dt = as.numeric(difftime(t,lastEventTime, units=unit))
    lastEventTime = t    
    
    #forgetting
    tn <- recomputeWeights(tn, dt, norm, tr)
    
    #add new link if not exist or set weight to 1 on existing link    
    if(tn[s,r]==0){
      tn <- tn + edge(s,r,weight=1.0, lt=0)
    }else{
      tn[s,r] <- 1
      E(tn,P=c(s,r))$lt <- 0
    }
    
    #count active nodes and components
    an_cn <- countActiveNodes(tn)
    #count weights and prepare histogram
    h <- countWeights(tn,bins)
    
    weightHist <- rbind(weightHist,h)
    anCnPerEvent <- rbind(anCnPerEvent,an_cn)
    trans <- c(trans, transitivity(tn))
    
  }
  
  results <- vector(mode = "list", length = 4)
  names(results) <- c("tn","an","wh","tr")
  results[[1]]=tn
  results[[2]]=anCnPerEvent
  results[[3]]=weightHist
  results[[4]]=trans
  return(results)
}


time = as.POSIXlt("2010-01-12 14:10:00")
unit = "days"
tr = 2
#select events from [time-threshold,time]
window <- tn_state(data,time,tr,unit)

#now make a nodes
tn <- findAllNodes(data)

#and create a network
window_results <- runEvents(tr = -tr, unit = unit, norm = 1, bins = seq(0,1,0.01), window = window, data = data)

subtn <- window_results[["tn"]]
plot.igraph(subtn, layout=layout.circle)