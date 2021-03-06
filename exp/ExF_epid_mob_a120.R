#parameters
dataset = "mobility"    #['mails','mobility'...]
norm = 1             #normalization of dt
unit = "secs"       #["secs","mins","hours","days","weeks","auto"] - do not use "auto" 
symdt = 0         #delay between ploting in second
show = FALSE          #plot or not
show_progress = FALSE 
threshold = exp(-120)     #level at which edge are deleted
bins = seq(0,1,0.01)
#bins = c(0,0.5,1)

library(igraph)
library(doParallel)
library(foreach)

registerDoParallel(cores=24)

loadData <- function(dataset){
  if(dataset=="mails"){
    data = read.csv2("e:/_Politechnika_Wroclawska/Projekty/ExF/manufacturing.csv")
    return(data)
  }
  if(dataset=="mobility"){
    data = read.csv2("/home/mkul/ms/mobility.csv")
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
  of = nrow(data)
  #for(i in 1:2000){
  for(i in 1:nrow(data)){      
    #read next event
    # print(paste(i,"of",of))
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


runInfection <- function(tr, norm, unit, symdt, show, bins,alpha,infect_start_step,gamma,critical_level){
  #pb = tkProgressBar("Events progress","%",1,nrow(data))
  #ib = tkProgressBar("Infection progress","%",0,1000)
  lastEventTime = as.POSIXct(data[1,3])
  
  tnet = tn
  
  niter = 10
  
  out_results = matrix(ncol=8)
  
  step = 0
  infected = 0
  infT = 0
  nth = 0
  start = TRUE
  
  tempV = V(tn)
  foreach(n=tempV)%dopar%{
    nth = nth + 1
    for(iter in 1:niter){
      tryCatch({
        #if(nth==1 && iter == 1) start = TRUE
        if(start){
          tn <- tnet
          inf_distr <- c()
          steps_time <- c()
          inf_time <- c()
          n_infected <- c()
          n_IS <- c()
          step = 0
          infected = 0
          for(i in 1:(nrow(data)-1)){
            step = step + 1
            ip=infectionProgress(tn)
            if(ip>critical_level){
              print("STOP!")
              break
            }
            if(step>limit){
              print("STOP! - out of limit")
              break
            }
            ip = round(ip*1000,digits=0)
            #read next event
            s <- toString(data[i,1])   #sender
            r <- toString(data[i,2])   #reciever
            t <- as.POSIXct(data[i,3]) #time (N_i)
            t_next <- as.POSIXct(data[i+1,3])
            steps_time <- c(steps_time,t)
            
            if(i==infect_start_step){
              #make nth node as bad one
              V(tn)[n]$infected = TRUE
              infected = infected + 1
              
              #progressbar
              #info <- paste(toString(ip)," infection done",sep="")
              #setTkProgressBar(ib,ip,label=info)
            }
            if(i>=infect_start_step){ #TODO:
              #randomly draw a time for next infection
              
              infT <- rexp(1,rate=gamma)
              t_0 <- t
              IS <- countIS(tn)
              inf_distr <- c(inf_distr,ip)
              n_IS <- c(n_IS,IS)
              
              time_to_next_inf = t_0 + (infT/IS)*alpha
              c_inf = 0
              while(time_to_next_inf<t_next){
                t_0 = time_to_next_inf
                #print(paste(time_to_next_inf,t_next,infT/IS*alpha))
                inf_time <- c(inf_time,time_to_next_inf)
                c_inf = c_inf + 1
                #infect random node:
                #take random infected node regardless of it has neighbors
                v = sample(V(tn)[V(tn)$infected==TRUE],1)
                if(length((V(tn)[neighbors(tn,v)])[V(tn)[neighbors(tn,v)]$infected==FALSE]$name)>0){
                  #and select ranom non infected neighbor
                  toinf = sample((V(tn)[neighbors(tn,v)])[V(tn)[neighbors(tn,v)]$infected==FALSE]$name,1)
                  V(tn)[toinf]$infected = TRUE
                  infected = infected + 1
                  
                  #progressbar
                  #info <- paste(toString(ip)," infection done",sep="")
                  #setTkProgressBar(ib,ip,label=info)
                }
                
                infT <- rexp(1,rate=gamma)
                IS <- countIS(tn)
                inf_distr <- c(inf_distr,ip)
                
                time_to_next_inf <- t_0 + infT/IS*alpha
              }
              n_infected <- c(n_infected,c_inf)
              delta = t_next - t_0
              infT = infT - delta*IS
              t_0 = t_next
              
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
            #info <- paste(toString(i),"/",toString(nrow(data))," done",sep="")
            #setTkProgressBar(pb,i, label=info)
            #show network
            if(show){
              plot.igraph(tn, vertex.size=5, layout=layout.sphere(tn), edge.width=E(tn)$weight*10)
            }
            Sys.sleep(symdt)
            
          }
          
          #close(pb)
          #close(ib)
          
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
          print(inf_distr)
          out_results <- rbind(out_results,results)
          path = paste("/home/grid/users/plgmkul/ms/res/1epid_full_sym_",toString(n),"_",toString(iter),"_",toString(step),sep="")
          save(out_results,file=path)
        }
      }, error = function(e){
        print(e)
        msg = paste(n,iter,sep=" ")
        loginfo(msg,logger="test")
      })
    }
  }
  return(out_results)
}


data = loadData(dataset)
tn = findAllNodes(data)

unit = "mins"
threshold = exp(-240)
alpha = 120
runInfection(threshold, norm, unit, symdt, show, bins, alpha, infect_start_step, gamma, critical_level)


#param = c(-260,-280,-300,-320,-340,-360,-380,-400,-420,-440,-460,-480,-500,-520,-540,-560,-580,-600)
#unit = "mins"
#
#tn = findAllNodes(data)
#
#foreach(p=param)%dopar%{
#  tn = findAllNodes(data)
#  threshold = exp(p)
#  res=runEvents(tr = threshold, norm = norm, unit = unit, symdt = symdt, show = show, bins = bins)
#  save(res,file = paste("/home/mkul/ms/res_mob/res_",p,"m",sep = ""))
#}

#an_cp <- results[["an"]]
#an <- an_cp[,1]
#cp <- an_cp[,2]
#dim(an) <- c(82928,1)
#dim(cp) <- c(82928,1)
#plot(an2,type='l',xlab='step',ylab='active nodes')
#plot(cp2,type='l',xlab='step',ylab='components')