#parameters
dataset = "mails"    #['mails',...]
norm = 1             #normalization of dt
unit = "hours"       #["secs","mins","hours","days","weeks","auto"] - do not use "auto" 
symdt = 0.1          #delay between ploting in second
show = FALSE          #plot or not
threshold = 0.01     #level at which 
bins = 100

loadData <- function(dataset){
  if(dataset=="mails"){
    data = read.csv2("manufacturing.csv")
    return(data)
  }
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

computeProbability <- function(){
  return(1)
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

runEvents <- function(tr, norm, unit, symdt, show, bins){
  pb = tkProgressBar("Progress","%",1,nrow(data))
  lastEventTime = as.POSIXct(data[1,3])
  
  anCnPerEvent = matrix(c(0,0),ncol=2)
  weightHist = matrix(ncol = bins)
  
  for(i in 1:nrow(data)){      
    #read next event    
    s = toString(data[i,1])
    r = toString(data[i,2])
    t = as.POSIXct(data[i,3])
    
    #time since last event
    dt = as.numeric(difftime(t,lastEventTime, units=unit))
    lastEventTime = t    
    
    #add new link if not exist or set weight to 1 on existing link    
    if(tn[s,r]==0){
      tn <- tn + edge(s,r,weight=1.0, lt=0)
    }else{
      tn[s,r] <- 0
      E(tn,P=c(s,r))$lt <- 0
    }
    
    #forgetting
    tn <- recomputeWeights(tn, dt, norm, tr)
    
    #count active nodes and components
    an_cn <- countActiveNodes(tn)
    #count weights and prepare histogram
    h <- countWeights(tn,bins)
    
    weightHist <- rbind(weightHist,h)
    anCnPerEvent <- rbind(anCnPerEvent,an_cn)
    
    #an <- temp[1]
    #cn <- temp[2]
    #print(paste("# of active nodes:",toString(an)," # of components:",toString(cn),sep=""))
    
    info <- paste(toString(i),"/",toString(nrow(data))," done",sep="")
    setTkProgressBar(pb,i, label=info)
    #show network
    if(show){
      plot.igraph(tn, vertex.size=5, layout=layout.sphere(tn), edge.width=E(tn)$weight*10)
    }
    Sys.sleep(symdt)
    
    #removing vertex - removes edges connected to removed vertex: tn <- delete.vertices(tn,'name')
    #removing edges - does not removing vertex: tn <- delete.edges(tn,E(tn,c("3","6")))
    
    
  }
  
  close(pb)
  return(c(tn,anCnPerEvent,weightHist))
}

#data <- loadData(dataset)
#tn <- findAllNodes(data)
results <- runEvents(threshold, norm, unit, symdt, show, bins)
