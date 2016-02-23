fib <- function(n){
  if (n == 0){
    return(0)
  }else{
    if (n == 1){
      return(1)
    }else{
      return (fib(n-1)+fib(n-2))
    }    
  }    
}

pb <- tkProgressBar("test progress bar", "Some information in %",
                    0, 100, 50)
Sys.sleep(0.5)
u <- c(0, sort(runif(20, 0, 100)), 100)
for(i in u) {
  Sys.sleep(0.1)
  info <- sprintf("%d%% done", round(i))
  setTkProgressBar(pb, i, sprintf("test (%s)", info), info)
}
Sys.sleep(5)
close(pb)

############################

max = 0
dts = matrix(c(0),ncol=1)
t = as.POSIXct(data[1,3])
pb <- tkProgressBar("Finding max difference time", "Progress",
                    0, nrow(data), 0)
for(i in 2:nrow(data)){
  lt = as.POSIXct(data[i,3])
  dt = as.numeric(difftime(lt,t, units='hours'))
  
  dts = rbind(dts,dt)
  
  t = lt
  info <- paste(toString(i),"/",toString(nrow(data))," done",sep="")
  setTkProgressBar(pb,i, label=info)
}
close(pb)

below = NULL
for(i in (1:length(an))){
  if(an[i]<50){
    below <- c(below,data[["EventDate"]][i])
  }
}



test1 <- function(a,b){
  c<-a+b
  return
}

test1(1,2)