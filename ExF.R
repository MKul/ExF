#parameters
dataset = "mails"    #['mails','mobility'...]
norm = 1             #normalization of dt
unit = "hours"       #["secs","mins","hours","days","weeks","auto"] - do not use "auto" 
symdt = 0         #delay between ploting in second
show = FALSE          #plot or not
show_progress = FALSE 
threshold = exp(-62.5)     #level at which edge are deleted
bins = seq(0,1,0.01)
#bins = c(0,0.5,1)


loadData <- function(dataset){
  if(dataset=="mails"){
    data = read.csv2("d:/_Politechnika_Wroclawska/Projekty/ExF/manufacturing.csv")
    return(data)
  }
  if(dataset=="mobility"){
    data = read.csv2("d:/_Politechnika_Wroclawska/Projekty/ExF/ExF/mobility.csv")
    data = data[,]
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
  if(show_progress){
    pb = tkProgressBar("Progress","%",1,nrow(data))
  }
  lastEventTime = as.POSIXct(data[1,3])
  
  anCnPerEvent = matrix(c(0,0),ncol=2)
  weightHist = matrix(ncol = length(bins))
  trans = c()
  
  #for(i in 1:2000){
  for(i in 1:nrow(data)){      
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
    
    if(show_progress){
      info <- paste(toString(i),"/",toString(nrow(data))," done",sep="")
      setTkProgressBar(pb,i, label=info)
    }
    
    #show network
    if(show){
      plot.igraph(tn, vertex.size=5, layout=layout.sphere(tn), edge.width=E(tn)$weight*10)
    }
    Sys.sleep(symdt)
    
  }
  if(show_progress){
    close(pb)
  }
  
  results <- vector(mode = "list", length = 4)
  names(results) <- c("tn","an","wh","tr")
  results[[1]]=tn
  results[[2]]=anCnPerEvent
  results[[3]]=weightHist
  results[[4]]=trans
  return(results)
}

#data <- loadData(dataset)
#tn <- findAllNodes(data)

#threshold = exp(-80)
#results1 <- runEvents(threshold, norm, unit, symdt, show, bins)
#threshold = exp(-168)
#results2 <- runEvents(threshold, norm, unit, symdt, show, bins)
#threshold = exp(-720)
#results3 <- runEvents(threshold, norm, unit, symdt, show, bins)

countIS <- function(tn){
  is = 0
  edges = get.edges(tn,E(tn))
  tryCatch({
    for(i in 1:ecount(tn)){
      v1 <- V(tn)[edges[i,1]]$infected
      v2 <- V(tn)[edges[i,2]]$infected
      print(toString(v1),toString(v2))
      if(xor(v1,v2)){
        is=is+1
      }
    }
    return(is)
  }, finally = {
    return(0)
  })
}

infectionProgress <- function(tn){
  n=0
  for(v in V(tn)){
    if(V(tn)[v]$infected==TRUE){
      n = n+1
    }
  }
  progress=(n/vcount(tn))
  return(progress)
}


runInfection <- function(tr, norm, unit, symdt, show, bins){
  if(show_progress){
    pb = tkProgressBar("Events progress","%",1,nrow(data))
    ib = tkProgressBar("Infection progress","%",1,vcount(tn))
  }
  lastEventTime = as.POSIXct(data[1,3])
  
  bad_one = 100     #event number when infection starts
  step = 0
  infected = 0
  gamma = 1
  infT = 0
  
  #for(i in 1:50){
  for(i in 1:(nrow(data)-1)){
    step = step + 1
    ip=infectionProgress(tn)
    if(ip>0.5){
      print("STOP!")
      break
    }
    #read next event    
    s = toString(data[i,1])   #sender
    r = toString(data[i,2])   #reciever
    t = as.POSIXct(data[i,3]) #time (N_i)
    t_next = as.POSIXct(data[i+1,3])
    
    if(i==bad_one){
      #make sender as bad one
      V(tn)[s]$infected = TRUE    
      infected = infected + 1
      if(show_progress){
        info <- paste(toString(infected),"/",toString(vcount(tn))," infection done",sep="")
        setTkProgressBar(ip,infected,label=info)
      }
    }
    if(i>=bad_one){ #TODO:
      #randomly draw a time for next infection
      infT = rexp(1,rate=gamma)
      t_0 = t
      IS = countIS(tn)
      
      if((t_0 + infT/IS)<t_next){
        t_0 = t_0+(infT/IS)
        #infect!
        for(v in V(tn)[V(tn)$infected==TRUE]){
          toinf = sample((V(tn)[neighbors(tn,v)])[V(tn)[neighbors(tn,v)]$infected==FALSE],1)
          V(tn)[toinf]$infected <- TRUE
          infected = infected + 1
          if(show_progress){
            info <- paste(toString(infected),"/",toString(vcount(tn))," done",sep="")
            setTkProgressBar(ip,infected,label=info)
          }
        }
        
      }else{
        delta = t_next - t_0
        infT = infT - delta*IS
        t_0 = t_next
      }
    }
    
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
    
    #output
    if(show_progress){
      info <- paste(toString(i),"/",toString(nrow(data))," done",sep="")
      setTkProgressBar(pb,i, label=info)
    }
    #show network
    if(show){
      plot.igraph(tn, vertex.size=5, layout=layout.sphere(tn), edge.width=E(tn)$weight*10)
    }
    Sys.sleep(symdt)
    
  }
  
  close(pb)
  results <- vector(mode = "list", length = 3)
  names(results) <- c("tn","an","wh")
  results[[1]]=tn
  results[[2]]=anCnPerEvent
  results[[3]]=weightHist
  return(results)
}

data = loadData(dataset)
tn = findAllNodes(data)
res2=runEvents(tr = threshold, norm = norm, unit = unit, symdt = symdt, show = show, bins = bins)
#runInfection(threshold, norm, unit, symdt, show, bins)

#an_cp <- results[["an"]]
#an <- an_cp[,1]
#cp <- an_cp[,2]
#dim(an) <- c(82928,1)
#dim(cp) <- c(82928,1)
#plot(an2,type='l',xlab='step',ylab='active nodes')
#plot(cp2,type='l',xlab='step',ylab='components')