#load data
dir = "d:/ExF/res_mob"
bresults = c()
for(file in list.files(path=dir)){
  print(paste(dir,file,sep="/"))
  load(paste(dir,file,sep="/"))
  bresults <- c(bresults,results)
}



#plot active nodes
for(i in 1:12){
  plot(bresults[[3*i-1]][,1],type="l")
}

#plot active components
for(i in 1:12){
  plot(bresults[[3*i-1]][,2],type="l")
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
