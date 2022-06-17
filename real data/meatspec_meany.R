#############################################################################
####-----------simulations for all methods-------------------------

rm(list = ls())
setwd(" ")##add your directory
source('loadAll.R')

####-------parameters for model---------------------------------------
num_simulations=20
data(meatspec)
N=nrow(meatspec)
P=ncol(meatspec)-1
X=meatspec[,1:P]
X=as.matrix(X)
y=meatspec[,(P+1)]
y=as.matrix(y)


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
}



