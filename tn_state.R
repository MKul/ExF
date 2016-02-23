tn_state <- function(d,ttime,threshold,unit){
  window <- c()
  begin <- F
  end <- F
  for(i in 1:length(d[,3])){
    t <- as.POSIXlt(d[i,3])
    diff <- as.numeric(difftime(ttime,t,units=unit))
    dexp <- exp(-diff)
    if(t<ttime && dexp>threshold){
      begin <- T
      window <- c(window,i)
    }else{
      if(begin) break
    }
  }
  window
}

make_network <- function(data, window){
  network <- graph.empty(directed=FALSE)
  for(i in 1:length(data[window,1])){
    s = toString(data[i,1])
    r = toString(data[i,2])
    #add new nodes if not exist
    if(length(which(V(network)$name == s))==0){
      network <- network + vertex(s, infected=FALSE)
      #print(paste("Node",s,"added.",sep=" "))
    }    
    if(length(which(V(network)$name == r))==0){
      network <- network + vertex(r, infected=FALSE)
      #print(paste("Node",r,"added.",sep=" "))
    }
    #add an edge if not exist
    if(length(E(network)[s %--% r])==0){
      network <- network + edge(s,r)
      #print(paste("Edge",s,r,"added.",sep=" "))
    }
    
  }
  network
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


recomputeWeights <- function(network, dt, norm, threshold){
  e <- ecount(network)
  if(e>0){
    todelete <- list()
    for(edge in E(network)){
      #new life time <- life time + normalized differece of event time 
      nlt <- E(network)[edge]$lt + (dt/norm)
      nw <- exp(-nlt)
      E(network)[edge]$weight <- nw
      tr <- threshold
      if(nw < tr){
        todelete <- append(todelete,edge)
      }else{
        E(network)[edge]$lt <- nlt
      }
    }  
    network <- delete.edges(network,todelete)
  }
  return(network)
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


runEvents <- function(tn, threshold, norm, unit, bins, window, data){
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
    difft = as.numeric(difftime(t,lastEventTime, units=unit))
    lastEventTime = t    
    
    #forgetting
    tn <- recomputeWeights(tn, difft, norm, threshold)
    
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


#time = as.POSIXlt("2010-01-12 14:10:00")
#unit = "days"
#tr = 2

#----select events from [time-threshold,time]----
#window <- tn_state(data,time,tr,unit)

#----now make a nodes----
#tn <- findAllNodes(data)

#----and create a network----
#window_results <- runEvents(tr = -tr, unit = unit, norm = 1, bins = seq(0,1,0.01), window = window, data = data)
#subtn <- window_results[["tn"]]
#plot.igraph(subtn, layout=layout.circle)

#----temporal clustering coefficients----
#begin <- as.POSIXct(data[1,"EventDate"])
#end <- as.POSIXct(data[length(data[,"EventDate"]),"EventDate"])
#bin_size <- as.numeric(end-begin)/100
#params <- seq(0,as.numeric(end-begin),as.numeric(end-begin)/100)
#wind_distr <- c()

#----events in windows----
#for(p in params){
#  p <- as.POSIXct(begin+(p*60*60*24))
#  print(p)
#  print(bin_size)
#  w <- tn_state(data = data, time = p, tr = bin_size, unit = "days")
#  wind_distr <- c(wind_distr,length(w))
#}

#----analysis----
#summary(wind_distr)
#hist(wind_distr)
#plot(sort(wind_distr),type="l")

#----clustering coefficient----
#clear_tn <- findAllNodes(data)
#globalcc <- c()
#localcc <- c()
#j <- 0
#for(p in params){
#  tn <- clear_tn
#  j <- j+1
#  print(j)
#  time <- as.POSIXct(begin+(p*60*60*24))
#  tr <- exp(-p*60*60*24)
#  w <- tn_state(data, time, bin_size, "days")
#  wr <- runEvents(tn, tr, 1, "days", seq(0,1,0.01), w, data)
  
#  subtn <- wr[["tn"]]
#  gcc <- transitivity(subtn,type="global")
#  lcc <- transitivity(subtn,type="local")
#  globalcc <- c(globalcc,gcc)
#  localcc <- c(localcc,list(lcc))
#}
#begin <- as.POSIXct(data[1,"t"])
#end <- as.POSIXct(data[length(data[,"t"]),"t"])
#bin_size <- as.numeric(end-begin)/100
#t <- as.POSIXct(data[200,"t"])
#u <- "mins"
#w <- tn_state(data, t, bin_size, u)

#t <- as.POSIXct(data[200,"t"])-as.difftime(exp(-240),units = "mins")

data <- loadData("mobility")
start <- data[1,3]
end <- data[length(data[,3]),3]
step <- difftime(end,start)/100
nets <- vector("list",100)
t <- as.POSIXlt(start)+step
for(i in 1:100){
  print(i)
  w <- tn_state(d = data, ttime = t,threshold = exp(-10), unit = "hours")
  t <- t + step
  network <- make_network(data = data,window = w)
  nets[[i]] <- network
}

# analysis

res <- read.csv2("d:/ExFMobilityNet/mobility100ExF_c.csv",sep=",")

exf <- as.vector(res[res[,2]==1,3])
deg <- as.vector(res[res[,2]==1,4])
ev <- as.vector(res[res[,2]==1,5])

all_exf <- as.vector(res[,3],mode="numeric")
all_deg <- as.vector(res[,4],mode="numeric")
all_ev <- as.vector(res[,5],mode="numeric")



jpeg(file="d:/exf_ev.jpg",width=800,height=600)
plot(all_exf,all_ev)
jpeg(file="d:/exf_deg.jpg",width=800,height=600)
plot(all_exf,all_deg)
jpeg(file="d:/deg_ev.jpg",width=800,height=600)
plot(all_deg,all_ev)
