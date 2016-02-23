#load data
dir = "d:/data/full_sym"
results = matrix(ncol=3)
for(file in list.files(path=dir)){
  print(paste(dir,file,sep="/"))
  load(paste(dir,file,sep="/"))
  for(i in nrow(out_results)){
    iter = out_results[i,c(1,2,3)]
    a = which(results[,1]==iter[1])
    b = which(results[,2]==iter[2])
    if(length(intersect(a,b))==0){
      results = rbind(results,do.call(cbind,iter))
    }
  }
  rm(out_results)
}



#plot everything
for(i in 2:length(out_results[,"inf_distr"])){
  plot(out_results[i,"inf_distr"][[1]]/(1000/167),type="l")
}

means = c()
for(i in 1:144){
  means = c(means,mean(results[results[,1]==i,3],na.rm=TRUE,trim = 0.1))
}

hist(means)

vars = c()
for(i in 1:144){
  vars = c(vars,var(results[results[,1]==i,3],na.rm=TRUE))
}

hist(vars,breaks=40)

sds = c()
for(i in 1:144){
  sds = c(sds,sd(results[results[,1]==i,3],na.rm=TRUE))
}

hist(sds)


epid_res_distr = cbind(1:144,means,sds,vars)
res_sort=epid_res_distr[sort.list(epid_res_distr[,2]),2:3]

twoord.plot(1:144,res_sort[,1],1:144,res_sort[,2],type="l",ylab="mean end step", rylab="standard deviation")
