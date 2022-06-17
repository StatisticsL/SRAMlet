#############################################################################
####-----------simulations for all methods-------------------------

rm(list = ls())
setwd("  ")####add your directory
source('loadAll.R')

####-------parameters for model---------------------------------------
num_simulations=100
Ntrain=2^10 #####a power of 2
n.train=floor(Ntrain*0.7)
n.val=Ntrain-n.train
Ntest=Ntrain
SD=1 ##standard deviation
SNR=3###the parameter in wave-functions
P=1000
S=4
if(S!=0){trueindex=seq(1:S);truenonindex=seq(1,P)[-trueindex]}
if(S==0){trueindex=NULL;truenonindex=seq(1,P)}

####----------parameters for sqrtAMlet,AMlet--------------------------
p0=2^0###a power of 2, the number of father wavelets
family="DaubExPhase"
filter.number=4
bc="periodic"
alpha=0.05
M_MC=1000
M_GEV=100
type='GEV'


# ####----------parameters for HAM--------------------------
# knots=floor(sqrt(n.train))
# ###for N=2^9,S=4,SNR=7,P=10,100,1000 and P=1000,N=2^10,
# lambda.pen  <- 0.955^(30:50)
# lambda.curv <- 2^c(-10:-18)

#####---------results save-------------------------------------
mse=matrix(NA,nrow = num_simulations,ncol = 8)
mse=as.data.frame(mse)
names(mse)=c('sqrtAMlet','oracleAMlet','madAMlet','HAM','SAM','mgcv','lasso1se','lassomin')
fdr.group.l0=mse
tpr.group.l0=mse

fitResults.sqrtAMlet=vector('list',length = num_simulations)
fitResults.oracleAMlet=vector('list',length = num_simulations)
fitResults.madAMlet=vector('list',length = num_simulations)
fitResults.HAM=vector('list',length = num_simulations)
fitResults.SAM=vector('list',length = num_simulations)
fitResults.mgcv=vector('list',length = num_simulations)
################################################################################
for (simu in 1:num_simulations){
  print(simu)
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
  
  x.train.part=x.train[1:n.train,]
  y.train.part=y.train[1:n.train]
  x.val.part=x.train[(n.train+1):Ntrain,]
  y.val.part=y.train[(n.train+1):Ntrain]
  
  # filename.data=paste('./simulation2/raw data/data','-N-',Ntrain,'-P-',P,'-S-',S,'-SNR-',SNR,'-simu-',simu,'.RData',sep = '')
  # save(x.train,y.train,x.train.part,y.train.part,x.val.part,y.val.part,x.test,y.test,mu.matrix.train,mu.matrix.test,
  #      file = filename.data)
  #############################################################################
  ####-------------all methods computation-----------------------------------
  
  ####----sqrtAMlet----------------------------------
  fit_sqrtAMlet=SRAMlet(y=y.train,X=x.train,lambda=NA,p0=p0,addHaar=F,addLinear=F,
                        filter.number=filter.number,family=family,
                        alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
  y.test.hat.sqrtAMlet=predict.sqrtAMlet(x.train,fit_sqrtAMlet,x.test)
  mse[simu,'sqrtAMlet']=mean((y.test-y.test.hat.sqrtAMlet$y.predict)^2)
  fdr.group.l0[simu,'sqrtAMlet']=groupL0Fdr(trueindex,truenonindex,fit_sqrtAMlet$index.selected)
  tpr.group.l0[simu,'sqrtAMlet']=groupL0Tpr(trueindex,truenonindex,fit_sqrtAMlet$index.selected)
  fitResults.sqrtAMlet[[simu]]=fit_sqrtAMlet
  
  # ####------oracle AMlet------------------------------------
  # fit_oracleAMlet=AMlet(y=y.train,X=x.train,lambda=NA,sigma=SD,p0=p0,filter.number=filter.number,family=family,bc=bc,
  #                      alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sqrtrss=T)
  # y.test.hat.oracleAMlet=predict.AMlet(x.train,fit_oracleAMlet,x.test)
  # mse[simu,'oracleAMlet']=mean((y.test-y.test.hat.oracleAMlet$y.predict)^2)
  # fdr.group.l0[simu,'oracleAMlet']=groupL0Fdr(trueindex,truenonindex,fit_oracleAMlet$index.selected)
  # tpr.group.l0[simu,'oracleAMlet']=groupL0Tpr(trueindex,truenonindex,fit_oracleAMlet$index.selected)
  # fitResults.oracleAMlet[[simu]]=fit_oracleAMlet
  
  ####-----mad AMlet-------------------------------------
  fit_madAMlet=AMlet.QUT.Mad(y=y.train,X=x.train,lambda=NA,sigma=NA,p0=p0,filter.number=filter.number,family=family,bc=bc,
                             alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
  y.test.hat.madAMlet=predict.AMlet(x.train,fit_madAMlet,x.test)
  mse[simu,'madAMlet']=mean((y.test-y.test.hat.madAMlet$y.predict)^2)
  fdr.group.l0[simu,'madAMlet']=groupL0Fdr(trueindex,truenonindex,fit_madAMlet$index.selected)
  tpr.group.l0[simu,'madAMlet']=groupL0Tpr(trueindex,truenonindex,fit_madAMlet$index.selected)
  fitResults.madAMlet[[simu]]=fit_madAMlet
  
  # ####----HAM--------------------------------------
  # fit_HAM=HAM(x.train=x.train.part,y.train=y.train.part,x.val=x.val.part,y.val=y.val.part,
  #             lambda.pen=lambda.pen,lambda.curv=lambda.curv,knots=knots, model = LinReg(),control = grpl.control(trace = 0))
  # mse[simu,'HAM']=fit_HAM$mse.Matrix[fit_HAM$index.pen,fit_HAM$index.curv]
  # fdr.group.l0[simu,'HAM']=groupL0Fdr(trueindex,truenonindex,fit_HAM$index.selected)
  # tpr.group.l0[simu,'HAM']=groupL0Tpr(trueindex,truenonindex,fit_HAM$index.selected)
  # fitResults.HAM[[simu]]=fit_HAM
  # 
  ####----SAM--------------------------------------
  fit_SAM=samQL(X=x.train.part,y=y.train.part)
  pred.SAM=predict(fit_SAM,newdata=x.val.part)
  lambda.SAM.index=which.min(colSums((pred.SAM$values-y.val.part)^2))
  index.selected.SAM=which(fit_SAM$func_norm[,lambda.SAM.index]!=0)
  pred.SAM=predict(fit_SAM,newdata = x.test)$values[,lambda.SAM.index]
  mse[simu,'SAM']=mean((pred.SAM-y.test)^2)
  fdr.group.l0[simu,'SAM']=groupL0Fdr(trueindex,truenonindex,index.selected.SAM)
  tpr.group.l0[simu,'SAM']=groupL0Tpr(trueindex,truenonindex,index.selected.SAM)
  fitResults.SAM[[simu]]=fit_SAM
  # 
  # ####----mgcv--------------------------------------
  # dat=as.data.frame(cbind(y.train,x.train))
  # names(dat)[1]='y.train'
  # fm <- paste('s(', names(dat)[-1], ')', sep = "", collapse = ' + ')
  # fm <- as.formula(paste('y.train ~', fm))
  # fit_mgcv=gam(fm,data=dat)
  # pred.mgcv=predict(fit_mgcv,newdata = as.data.frame(cbind(y.test,x.test)),type = 'terms')
  # y.test.pred.mgcv=rowSums(pred.mgcv)
  # mse[simu,'mgcv']=mean((y.test-y.test.pred.mgcv)^2)
  # ll=length(fit_mgcv$coefficients[-1])
  # coeff.matrix=matrix(fit_mgcv$coefficients[-1],nrow=ll/P,ncol = P)
  # index.selected.mgcv=which(colSums(abs(coeff.matrix))!=0)
  # fdr.group.l0[simu,'mgcv']=groupL0Fdr(trueindex,truenonindex,index.selected.mgcv)
  # tpr.group.l0[simu,'mgcv']=groupL0Tpr(trueindex,truenonindex,index.selected.mgcv)
  # fitResults.mgcv[[simu]]=fit_mgcv

  lambda=cv.glmnet(x.train,y.train)
  fit=glmnet(x.train,y.train,family = 'gaussian',lambda = lambda$lambda.min)
  fdr.group.l0[simu,'lassomin']=groupL0Fdr(trueindex,truenonindex,which(fit$beta!=0))
  tpr.group.l0[simu,'lassomin']=groupL0Tpr(trueindex,truenonindex,which(fit$beta!=0))
  y.test.pre.min=predict(lambda,x.test,s='lambda.min')
  mse[simu,'lassomin']=mean((y.test-y.test.pre.min)^2)
  ########################################################
  fit=glmnet(x.train,y.train,family = 'gaussian',lambda = lambda$lambda.1se)
  fdr.group.l0[simu,'lasso1se']=groupL0Fdr(trueindex,truenonindex,which(fit$beta!=0))
  tpr.group.l0[simu,'lasso1se']=groupL0Tpr(trueindex,truenonindex,which(fit$beta!=0))
  y.test.pre.1se=predict(lambda,x.test,s='lambda.1se')
  mse[simu,'lasso1se']=mean((y.test-y.test.pre.1se)^2)
}

results=matrix(NA,nrow = 8,ncol = 6)
results=as.data.frame(results)
row.names(results)=c('sqrtAMlet','oracleAMlet','madAMlet','HAM','SAM','mgcv','lassomin','lasso1se')
names(results)=c('mse.mean','mse.sd','groupl0fdr.mean','groupl0fdr.sd','groupl0tpr.mean','groupl0tpr.sd')
results[,'mse.mean']=apply(mse,2,mean)
results[,'mse.sd']=apply(mse,2,sd)/sqrt(num_simulations)
results[,'groupl0fdr.mean']=apply(fdr.group.l0,2,mean)
results[,'groupl0fdr.sd']=apply(fdr.group.l0,2,sd)/sqrt(num_simulations)
results[,'groupl0tpr.mean']=apply(tpr.group.l0,2,mean)
results[,'groupl0tpr.sd']=apply(tpr.group.l0,2,sd)/sqrt(num_simulations)

save(mse,fdr.group.l0,tpr.group.l0,results,file = paste('./simulation2/results/compution-N-',Ntrain,'-P-',P,'-S-',S,'-SNR-',SNR,'.RData',sep = ''))
save(fitResults.sqrtAMlet,fitResults.oracleAMlet,fitResults.madAMlet,fitResults.HAM,fitResults.SAM,fitResults.mgcv,
     file =paste( './simulation2/results/fitting-N-',Ntrain,'-P-',P,'-S-',S,'-SNR-',SNR,'.RData',sep = ''))





