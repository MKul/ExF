#parameters
dataset = "mails"    #['mails','mobility'...]
norm = 1             #normalization of dt
unit = "hours"       #["secs","mins","hours","days","weeks","auto"] - do not use "auto" 
symdt = 0         #delay between ploting in second
show = FALSE          #plot or not
show_progress = FALSE 
threshold = exp(-720)     #level at which edge are deleted
bins = seq(0,1,0.01)
infect_start_step = 1000  #event number when infection starts
critical_level = 0.4
alpha = 60*60*24
gamma = 1
#limit = 20000
res_files_name = "epid_full_sym_man_"
args <- commandArgs(trailingOnly = TRUE)
node <- strtoi(args[1])
#node <- 164
res_dir <- "/home/mkul/ms/res_man_epid_a30/"
#res_dir <- "d:/epid_sym_man_tmp/"

library(igraph)
library(doParallel)

registerDoParallel(cores=1)
#registerDoParallel(cores=8)
print("Running...")

loadData <- function(dataset){
  if(dataset=="mails"){
    data = read.csv2("/home/mkul/ms/manufacturing.csv")
    #data = read.csv2("d:/_Politechnika_Wroclawska/Projekty/ExF/manufacturing.csv")
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
    }    
    if(length(which(V(tn)$name == r))==0){
      tn <- tn + vertex(r, infected=FALSE)
    }
  }
  adjMat <- matrix(-10,nrow = length(V(tn)), ncol = length(V(tn)))
  tn <- permute.vertices(tn,V(tn)$name)
  return(list(tn,adjMat))
}


recomputeWeights <- function(n, dt, norm, tr){
  adjMat <- adjMat * exp(-dt)
  v1 <- which(adjMat < tr)
  v2 <- which(adjMat > 0)
  todel <- intersect(v1,v2)
  
  ids = c()
  todelete = c()
  dimension = dim(adjMat)[1]
  
  #for(td in todel){  
  foreach(td=todel)%dopar%{
    x = td%%dimension
    if(x==0){
      x=dimension
      y = td%/%dimension
    }else{
      y = td%/%dimension+1
    }
    ids = c(ids, get.edge.ids(n,c(x,y)))
    adjMat[td] = -10
  }  
  #print(E(n))
  ids = unique(ids[which(ids!=0)])
  #print(ids)
  n <- delete.edges(n,ids)
  
  return(list(n,adjMat))
}

countWeights <- function(n,bins){
  if(ecount(n)>0){
    h <- hist(E(n)$weight,breaks=bins,plot=FALSE)
    return(h$counts)
  }else{
    return(matrix(0,ncol=bins))
  }
}

countActiveNodes <- function (n){
  an <- 0
  components <- decompose.graph(n)
  cn <- 0
  for(x in components){
    if(ecount(x)>0){
      an <- an + vcount(x)
      cn <- cn + 1
    }
  }
  return(c(an,cn))
}


runEvents <- function(tr, norm, unit, symdt, show, bins){
  
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


countIS <- function(tn){
  is <- 0
  edges <- get.edges(tn,E(tn))
  if(length(edges)>0){
    for(i in 1:ecount(tn)){
      v1 <- V(tn)[edges[i,1]]$infected
      v2 <- V(tn)[edges[i,2]]$infected
      if(xor(v1,v2)){
        is=is+1
      }
    }
  }
  return(is)
}

infectionProgress <- function(tn){
  n <- 0
  for(v in V(tn)){
    if(V(tn)[v]$infected==TRUE){
      n <- n+1
    }
  }
  progress <- (n/vcount(tn))
  return(progress)
}


runInfection <- function(tr, norm, unit, symdt, show, bins,alpha,infect_start_step,gamma,critical_level, tn, adjMat){
  
  lastEventTime = as.POSIXct(data[1,3])
  tnet = tn
  aMat = adjMat
  niter = 10
  results_man = matrix(ncol=8)  
  infT = 0
  tempV = V(tn)
  number_of_nodes = vcount(tn)
  
  #foreach(n=tempV)%dopar%{
  n = node
  
  for(iter in 1:niter){
    start = TRUE
    if(start){
      tn <- tnet
      adjMat <- aMat
      inf_distr <- c()
      steps_time <- c()
      inf_time <- c()
      n_infected <- c()
      n_IS <- c()
      step = 0
      infected = 0
      inf_started = FALSE
      ip = 0
      stime=Sys.time()
      # iterate through the events
      for(i in 1:(nrow(data)-1)){
        tryCatch({
          if(i%%20 == 0 && inf_started){
            print(paste(i,"progress:",ip))
            print(Sys.time()-stime)
            stime=Sys.time()
            #V(tn)$color = V(tn)$infected
            #V(tn)$color = gsub("TRUE", "green", V(tn)$color)
            #V(tn)$color = gsub("FALSE", "gray", V(tn)$color)
            #plot.igraph(tn, vertex.size = 5, layout = layout.circle)
          }
          
          step = step + 1
          if(inf_started){
            #ip=infectionProgress(tn)
            ip = infected/number_of_nodes
            #print(paste(i,"Progress:",ip,"of",critical_level,sep=" "))
            if(ip>critical_level){
              print("STOP!")
              break
            }
            #if(step>limit){
            #  print("STOP! - out of limit")
            #  break
            #}
          }
          #read next event
          sender <- toString(data[i,1])   #sender
          reciever <- toString(data[i,2])   #reciever
          t <- as.POSIXct(data[i,3]) #time (N_i)
          t_next <- as.POSIXct(data[i+1,3])
          steps_time <- c(steps_time,t)
          
          if(i==infect_start_step){
            inf_started = TRUE
            V(tn)[n]$infected = TRUE
            infected = infected + 1
            print("Infection began.")
          }
          if(i>=infect_start_step){
            #randomly draw a time for next infection
            tempt = Sys.time()
            infT <- rexp(1,rate=gamma)
            t_0 <- t
            IS <- countIS(tn)
            inf_distr <- c(inf_distr,ip)
            n_IS <- c(n_IS,IS)
            
            time_to_next_inf = t_0 + (infT/IS)*alpha
            c_inf = 0
            tryCatch({
              while(time_to_next_inf<t_next && IS>0){
                t_0 = time_to_next_inf
                inf_time <- c(inf_time,time_to_next_inf)
                c_inf = c_inf + 1
                #infect random node:
                #take random infected node regardless of it has neighbors
                v = sample(V(tn)[V(tn)$infected==TRUE],1)
                #print("Im trying")
                #print(V(tn)[neighbors(tn,v)]$infected)
                if(length((V(tn)[neighbors(tn,v)])[V(tn)[neighbors(tn,v)]$infected==FALSE])>0){
                  #and select random non infected neighbor
                  toinf = sample((V(tn)[neighbors(tn,v)])[V(tn)[neighbors(tn,v)]$infected==FALSE],1)
                  V(tn)[toinf]$infected = TRUE
                  infected = infected + 1
                  print('infect!')
                }
                
                infT <- rexp(1,rate=gamma)
                IS <- countIS(tn)
                inf_distr <- c(inf_distr,ip)
                
                time_to_next_inf <- t_0 + infT/IS*alpha
              }
            },error = function(e){
              print(paste("No such vertex:",v))
            })
            n_infected <- c(n_infected,c_inf)
            delta = t_next - t_0
            infT = infT - delta*IS
            t_0 = t_next
            
          }
          #time since last event
          dt = as.numeric(difftime(t,lastEventTime, units=unit))  
          lastEventTime = t         
          
          #forgetting
          if(dt > 0){
            out <- recomputeWeights(tn, dt, norm, tr)
            tn <- out[[1]]
            adjMat <- out[[2]]
          }
          
          #add new link if not exist or set weight to 1 on existing link
          sender = as.integer(sender)
          reciever = as.integer(reciever)
          if(tn[sender,reciever]==0){
            tn <- tn + edge(sender,reciever)
            adjMat[sender,reciever] <- 1.0
            adjMat[reciever,sender] <- 1.0
          }else{
            adjMat[sender,reciever] <- 1.0
            adjMat[reciever,sender] <- 1.0
          }
        },error = function(e){
          print(paste("Error on event:",i,"s:",sender,"r",reciever))
          print(e)
        })
      }
      #saving results
      results <- vector(mode = "list", length = 8)
      names(results) <- c("node","iter","end_step","inf_distr","steps_time","inf_time","number_of_infected","number_of_IS")
      results[[1]]=n
      results[[2]]=iter
      results[[3]]=step
      results[[4]]=inf_distr
      results[[5]]=steps_time
      results[[6]]=inf_time
      results[[7]]=n_infected
      results[[8]]=n_IS
      
      results_man <- rbind(results_man,results)
      path = paste(res_dir,res_files_name,toString(n),"_",toString(iter),"_",toString(step),sep="")
      save(results_man,file=path)
      results_man = matrix(ncol=8)
    }
  }
}

data <- loadData(dataset)
out <- findAllNodes(data)
tn <- out[[1]]
adjMat = out[[2]]

print(paste("Infect from node ",node))
tryCatch({
  runInfection(threshold, norm, unit, symdt, show, bins, alpha, infect_start_step, gamma, critical_level, tn, adjMat = adjMat)
}, error = function(e){
  print(e)
  traceback()
})
