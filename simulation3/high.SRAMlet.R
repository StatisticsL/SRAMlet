#############################################################################
####-----------simulations for all methods-------------------------

rm(list = ls())
t1=Sys.time()
setwd("  ")##add your directory
source('loadAll.R')

####-------parameters for model---------------------------------------
j=12
num_simulations=1
Ntrain=2^(j) #####a power of 2
P=2^(j+1)
# P=10
Ntest=Ntrain
SD=1 ##standard deviation
SNR=3###the parameter in wave-functions
S=4
if(S!=0){trueindex=seq(1:S);truenonindex=seq(1,P)[-trueindex]}
if(S==0){trueindex=NULL;truenonindex=seq(1,P)}

####----------parameters for sqrtAMlet,AMlet--------------------------
p0=2^0###a power of 2, the number of father wavelets
J=log2(p0)
family="DaubExPhase"
filter.number=4
bc="periodic"
alpha=0.05
M_MC=1000
M_GEV=100
type='GEV'

#####---------results save-------------------------------------
mse=matrix(NA,nrow = num_simulations,ncol = 2)
mse=as.data.frame(mse)
names(mse)=c('sqrtAMlet','madAMlet')
fdr.group.l0=mse
tpr.group.l0=mse
sigma.all=mse

fitResults.sqrtAMlet=vector('list',length = num_simulations)
fitResults.oracleAMlet=vector('list',length = num_simulations)
fitResults.madAMlet=vector('list',length = num_simulations)

################################################################################
for (simu in 1:num_simulations){
  ####-------------generate data-------------------------------------------
  x.train=matrix(runif(Ntrain*P),Ntrain,P)
  x.test=matrix(runif(Ntest*P),Ntest,P)
  
  mu.matrix.train=matrix(0,nrow = Ntrain,ncol = S)
  mu.matrix.test=matrix(0,nrow = Ntest,ncol=S)
  
  mu.matrix.train[,1]=waveblocks(x=x.train[,1],snr = SNR)
  mu.matrix.test[,1]=waveblocks(x=x.test[,1],snr = SNR)
  mu.matrix.train[,2]=wavedoppler(x=x.train[,2],snr = SNR)
  mu.matrix.test[,2]=wavedoppler(x=x.test[,2],snr = SNR)
  mu.matrix.train[,3]=wavebumps(x=x.train[,3],snr = SNR)
  mu.matrix.test[,3]=wavebumps(x=x.test[,3],snr = SNR)
  mu.matrix.train[,4]=waveheavisine(x=x.train[,4],snr = SNR)
  mu.matrix.test[,4]=waveheavisine(x=x.test[,4],snr = SNR)
  
  y.train=apply(mu.matrix.train,1,sum)+rnorm(Ntrain,sd=SD)
  y.test=apply(mu.matrix.test,1,sum)
  
  # filename.data=paste('./simulation3/raw data/data','-N-',Ntrain,'-P-',P,'-S-',S,'-SNR-',SNR,'-simu-',simu,'.RData',sep = '')
  # save(x.train,y.train,x.test,y.test,mu.matrix.train,mu.matrix.test,
  #      file = filename.data)
  ####----sqrtAMlet----------------------------------
  fit_sqrtAMlet=SRAMlet(y=y.train,X=x.train,lambda=NA,p0=p0,addHaar = F,addLinear = F,
                        filter.number=filter.number,family=family,
                        alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sqrtrss=T)
  y.test.hat.sqrtAMlet=predict.sqrtAMlet(x.train,fit_sqrtAMlet,x.test)
  mse[simu,'sqrtAMlet']=mean((y.test-y.test.hat.sqrtAMlet$y.predict)^2)
  fdr.group.l0[simu,'sqrtAMlet']=groupL0Fdr(trueindex,truenonindex,fit_sqrtAMlet$index.selected)
  tpr.group.l0[simu,'sqrtAMlet']=groupL0Tpr(trueindex,truenonindex,fit_sqrtAMlet$index.selected)
  fitResults.sqrtAMlet[[simu]]=fit_sqrtAMlet
  sigma.all[simu,'sqrtAMlet']=mean(fit_sqrtAMlet$sqrtrss.vec)

  # # 
  # ####-----mad AMlet-------------------------------------
  # fit_madAMlet=AMlet.QUT.Mad(y=y.train,X=x.train,lambda=NA,sigma=NA,p0=p0,filter.number=filter.number,family=family,bc=bc,
  #                            alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sqrtrss=T)
  # y.test.hat.madAMlet=predict.AMlet(x.train,fit_madAMlet,x.test)
  # mse[simu,'madAMlet']=mean((y.test-y.test.hat.madAMlet$y.predict)^2)
  # fdr.group.l0[simu,'madAMlet']=groupL0Fdr(trueindex,truenonindex,fit_madAMlet$index.selected)
  # tpr.group.l0[simu,'madAMlet']=groupL0Tpr(trueindex,truenonindex,fit_madAMlet$index.selected)
  # fitResults.madAMlet[[simu]]=fit_madAMlet
  # sigma.all[simu,'madAMlet']=fit_madAMlet$sigma
  
  print(paste('Completed',simu,sep = ' '))
}

results=matrix(NA,nrow = 2,ncol = 6)
results=as.data.frame(results)
row.names(results)=c('sqrtAMlet','madAMlet')
names(results)=c('mse.mean','mse.sd','groupl0fdr.mean','groupl0fdr.sd','groupl0tpr.mean','groupl0tpr.sd')
results[,'mse.mean']=apply(mse,2,mean)
results[,'mse.sd']=apply(mse,2,sd)
results[,'groupl0fdr.mean']=apply(fdr.group.l0,2,mean)
results[,'groupl0fdr.sd']=apply(fdr.group.l0,2,sd)
results[,'groupl0tpr.mean']=apply(tpr.group.l0,2,mean)
results[,'groupl0tpr.sd']=apply(tpr.group.l0,2,sd)

save(mse,fdr.group.l0,tpr.group.l0,results,sigma.all,file = paste('./simulation3/results/compution-N-',Ntrain,'-P-',P,'-S-',S,'-SNR-',SNR,'.RData',sep = ''))
save(fitResults.sqrtAMlet,fitResults.madAMlet,
     file =paste( './simulation3/results/fitting-N-',Ntrain,'-P-',P,'-S-',S,'-SNR-',SNR,'.RData',sep = ''))
t2=Sys.time()

print(paste('N=',Ntrain,'P=',P,sep = ' '))
print(results)
print(t2-t1)


