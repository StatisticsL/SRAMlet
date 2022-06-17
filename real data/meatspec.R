#############################################################################
####-----------simulations for all methods-------------------------

rm(list = ls())
setwd(" ")##add your directory
source('loadAll.R')

####-------parameters for model---------------------------------------
num_simulations=1
data(meatspec)
N=nrow(meatspec)
P=ncol(meatspec)-1
X=meatspec[,1:P]
X=as.matrix(X)
y=meatspec[,(P+1)]
y=as.matrix(y)

###not scale
#X=scale(X)
# l2norm=sqrt(apply(X^2,2,sum))
# X=X/matrix(rep(l2norm,N),N,P,byrow = TRUE)


####----------parameters for sqrtAMlet,AMlet--------------------------
p0=2^0###a power of 2, the number of father wavelets
family="DaubExPhase"
filter.number=3
bc="periodic"
alpha=0.05
M_MC=1000
M_GEV=100
type='GEV'


####----------parameters for HAM--------------------------
knots=floor(sqrt(89))
###for N=2^9,S=4,SNR=7,P=10,100,1000 and P=1000,N=2^10,
# lambda.pen  <- 0.955^(30:50)
lambda.pen  <- 0.955^(30:50)
lambda.curv <- 2^c(-10:-18)

#####---------results save-------------------------------------
mse=matrix(NA,nrow = num_simulations,ncol = 5)
mse=as.data.frame(mse)
names(mse)=c('sqrtAMlet','madAMlet','HAM','SAM','mgcv')

modelsize=mse
selectedindex.SRAMlet=vector('list',length = num_simulations)
selectedindex.madAMlet=vector('list',length = num_simulations)
selectedindex.HAM=vector('list',length = num_simulations)
selectedindex.SAM=vector('list',length = num_simulations)
selectedindex.mgcv=vector('list',length = num_simulations)

fitResults.sqrtAMlet=vector('list',length = num_simulations)
fitResults.madAMlet=vector('list',length = num_simulations)
fitResults.HAM=vector('list',length = num_simulations)
fitResults.SAM=vector('list',length = num_simulations)
fitResults.mgcv=vector('list',length = num_simulations)
################################################################################
for (simu in 1:num_simulations){
  train.index=sample(seq(1,N),size = 2^7,replace = FALSE)
  test.index=setdiff(seq(1,N),train.index)
  
  x.train=X[train.index,]
  x.test=X[test.index,]
  y.train=y[train.index]
  y.test=y[test.index]
 
  train.part.index=sample(train.index,size =89 ,replace = FALSE)
  val.part.index=setdiff(train.index,train.part.index)
  
  x.train.part=X[train.part.index,]
  y.train.part=y[train.part.index]
  x.val.part=X[val.part.index,]
  y.val.part=y[val.part.index]
  #############################################################################
  ####-------------all methods computation-----------------------------------
  
  ####----sqrtAMlet----------------------------------
  fit_sqrtAMlet=SRAMlet(y=y.train,X=x.train,lambda=NA,p0=p0,addHaar=T,addLinear=T,filter.number=filter.number,family=family,
                          alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
  y.test.hat.sqrtAMlet=predict.sqrtAMlet(x.train,fit_sqrtAMlet,x.test)
  mse[simu,'sqrtAMlet']=mean((y.test-y.test.hat.sqrtAMlet$y.predict)^2)
  selectedindex.SRAMlet[[simu]]=fit_sqrtAMlet$index.selected
  modelsize[simu,'sqrtAMlet']=length(fit_sqrtAMlet$index.selected)
  fitResults.sqrtAMlet[[simu]]=fit_sqrtAMlet
  
  ####-----mad AMlet-------------------------------------
  fit_madAMlet=AMlet(y=y.train,X=x.train,lambda=NA,sigma=NA,p0=p0,addHaar=T,addLinear=F,filter.number=filter.number,family=family,
                             alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
  y.test.hat.madAMlet=predict.AMlet(x.train,fit_madAMlet,x.test)
  mse[simu,'madAMlet']=mean((y.test-y.test.hat.madAMlet$y.predict)^2)
  modelsize[simu,'madAMlet']=length(fit_madAMlet$index.selected)
  selectedindex.madAMlet[[simu]]=fit_madAMlet$index.selected
  fitResults.madAMlet[[simu]]=fit_madAMlet
  #
  # 
  ####----HAM--------------------------------------
  fit_HAM=HAM(x.train=x.train.part,y.train=y.train.part,x.val=x.val.part,y.val=y.val.part,
              lambda.pen=lambda.pen,lambda.curv=lambda.curv,knots=knots, model = LinReg(),control = grpl.control(trace = 0))
  mse[simu,'HAM']=fit_HAM$mse.Matrix[fit_HAM$index.pen,fit_HAM$index.curv]
  modelsize[simu,'HAM']=length(fit_HAM$index.selected)
  selectedindex.HAM[[simu]]=fit_HAM$index.selected
  fitResults.HAM[[simu]]=fit_HAM

  ####----SAM------------------------------###problem!!!
  fit_SAM=samQL(X=x.train.part,y=y.train.part)
  pred.SAM=predict(fit_SAM,newdata=x.val.part)
  lambda.SAM.index=which.min(colSums((pred.SAM$values-y.val.part)^2))
  index.selected.SAM=which(fit_SAM$func_norm[,lambda.SAM.index]!=0)
  pred.SAM=predict(fit_SAM,newdata = x.test)$values[,lambda.SAM.index]
  mse[simu,'SAM']=mean((pred.SAM-y.test)^2)
  modelsize[simu,'SAM']=length(index.selected.SAM)
  selectedindex.SAM[[simu]]=index.selected.SAM
  fitResults.SAM[[simu]]=fit_SAM
  
  # ####----mgcv----------problem: the dimension p is too large-------------------------
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
  # modelsize[simu,'mgcv']=length( index.selected.mgcv )
  # selectedindex.mgcv[[simu]]=index.selected.mgcv
  # fitResults.mgcv[[simu]]=fit_mgcv
  
  print(paste('Completed ',simu,sep = ''))
}


results=matrix(NA,nrow = 5,ncol = 4)
results=as.data.frame(results)
row.names(results)=c('sqrtAMlet','madAMlet','HAM','SAM','mgcv')
names(results)=c('mse.mean','mse.sd','size.mean','size.sd')
results[,'mse.mean']=apply(mse,2,mean)
results[,'mse.sd']=apply(mse,2,sd)
results[,'size.mean']=apply(modelsize,2, mean)
results[,'size.sd']=apply(modelsize,2, sd)
# 
# save(mse,modelsize,results,selectedindex.SRAMlet,selectedindex.madAMlet,selectedindex.HAM,selectedindex.SAM
#      ,selectedindex.mgcv,file = './simulation_realdata/results-meatspec.RData')
# save(fitResults.sqrtAMlet,fitResults.madAMlet,fitResults.HAM,fitResults.SAM,fitResults.mgcv,
#      file ='./simulation_realdata/fitting-meatspec.RData')
# 
# 


