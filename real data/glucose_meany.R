#############################################################################
####-----------simulations for all methods-------------------------

rm(list = ls())
setwd(" ")##add your directory
source('loadAll.R')
####-------parameters for model---------------------------------------
num_simulations=20
data(NIR)
X=NIR$xNIR
y=NIR$yGlcEtOH$Glucose
X=as.matrix(X)
y=as.matrix(y)
N=nrow(X)##N=166  training sets:128 test sets:38
P=ncol(X)## P=235

#####---------results save-------------------------------------
mse=matrix(NA,nrow = num_simulations,ncol = 1)
mse=as.data.frame(mse)
names(mse)='meany'

################################################################################
for (simu in 1:num_simulations){
  train.index=sample(seq(1,N),size = 2^7,replace = FALSE)
  test.index=setdiff(seq(1,N),train.index)
  
  x.train=X[train.index,]
  x.test=X[test.index,]
  y.train=y[train.index]
  y.test=y[test.index]
  
  meany=rep(mean(y.train),length(y.test))
  mse[simu,'meany']=mean((meany-y.test)^2)
  #############################################################################
}

mean(mse[,'meany'])
sd(mse[,'meany'])



