#############################################################################
####-----------simulations for all methods-------------------------

rm(list = ls())
setwd(" ")##add your directory
source('loadAll.R')
####-------parameters for model---------------------------------------
num_simulations=20
data(NIR)
X=NIR$xNIR
y=NIR$yGlcEtOH$Ethanol
X=as.matrix(X)
y=as.matrix(y)
N=nrow(X)##N=166  training sets:128 test sets:38
P=ncol(X)## P=235

#####---------results save-------------------------------------
mse=matrix(NA,nrow = num_simulations,ncol = 2)
mse=as.data.frame(mse)
names(mse)=c('lassomin','lasso1se')
modelsize=mse
selectedindex.min=vector('list',length = num_simulations)
selectedindex.1se=vector('list',length = num_simulations)


################################################################################
for (simu in 1:num_simulations){
  train.index=sample(seq(1,N),size = 2^7,replace = FALSE)
  test.index=setdiff(seq(1,N),train.index)
  
  x.train=X[train.index,]
  x.test=X[test.index,]
  y.train=y[train.index]
  y.test=y[test.index]
  
  lambda=cv.glmnet(x.train,y.train)
  fit=glmnet(x.train,y.train,family = 'gaussian',lambda = lambda$lambda.min)
  selectedindex.min[[simu]]=which(fit$beta!=0)
  modelsize[simu,'lassomin']=length(which(fit$beta!=0))
  y.test.pre.min=predict(lambda,x.test,s='lambda.min')
  mse[simu,'lassomin']=mean((y.test-y.test.pre.min)^2)
  ########################################################
  fit=glmnet(x.train,y.train,family = 'gaussian',lambda = lambda$lambda.1se)
  selectedindex.1se[[simu]]=which(fit$beta!=0)
  modelsize[simu,'lasso1se']=length(which(fit$beta!=0))
  y.test.pre.1se=predict(lambda,x.test,s='lambda.1se')
  mse[simu,'lasso1se']=mean((y.test-y.test.pre.1se)^2)
  #############################################################################
}

apply(mse, 2, mean)
apply(mse, 2, sd)
apply(modelsize, 2, mean)
apply(modelsize, 2, sd)





