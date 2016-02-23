#load data
dir = "d:/results_R/ExF/epid_sym_mob"
bresults = c()

ends <- vector("list",274)
for(file in list.files(path=dir)){
  load(paste(dir,file,sep="/"))
  node <- strtoi(substr(file,gregexpr(pattern = "_", file)[[1]][4]+1,gregexpr(pattern = "_", file)[[1]][5]-1))
  end <- strtoi(substr(file,gregexpr(pattern = "_", file)[[1]][6]+1,nchar(file)))
  ends[[node]] <- c(ends[[node]],end)
  bresults <- c(bresults,results_mob_a30[2,])
}

# plot outcomes
means = c()
vars = c()
sds = c()
for(n in ends){
  if(!is.null(n)){
    means <- c(means, mean(n))
    vars <- c(vars, var(n))
    sds <- c(sds, sd(n))
  }else{
    means <- c(means, "NA")
    vars <- c(vars, "NA")
    sds <- c(sds, "NA")
  }
}

# ExF
res <- read.csv2("d:/results_R/ExF/ExFMobilityNet/mobility100ExF_c.csv",sep=",")

all_exf <- as.vector(res[,3],mode="numeric")
all_deg <- as.vector(res[,4],mode="numeric")
all_ev <- as.vector(res[,5],mode="numeric")
exfs <- c()
degs <- c()
evs <- c()
for(i in 1:274){
  exf <- mean(as.vector(res[res[,2]==i,3],mode="numeric"),na.rm = TRUE)
  deg <- mean(as.vector(res[res[,2]==i,4],mode="numeric"),na.rm = TRUE)
  ev <- mean(as.vector(res[res[,2]==i,5],mode="numeric"),na.rm = TRUE)
  exfs <- c(exfs,exf)
  degs <- c(degs,deg)
  evs <- c(evs,ev)
}

epid_res_distr = cbind(1:length(means),as.numeric(means),as.numeric(sds),as.numeric(vars),as.numeric(exfs),as.numeric(degs),as.numeric(evs))
res_sort=epid_res_distr[sort.list(epid_res_distr[,2]),1:7]


#twoord.plot(1:length(means),res_sort[,2],1:length(means),res_sort[,5],
#            main="End steps per node vs ExF", type="l",ylab="mean", rylab="exf",ylog = TRUE)
#twoord.plot(1:length(means),res_sort[,2],1:length(means),res_sort[,6],
#            main="End steps per node vs degree", type="l",ylab="mean", rylab="deg",ylog = TRUE)
#twoord.plot(1:length(means),res_sort[,2],1:length(means),res_sort[,7],
#            main="End steps per node vs eigenvalue", type="l",ylab="mean", rylab="ev",ylog = TRUE)
#axis(side=1,at = seq(0,300,10))


