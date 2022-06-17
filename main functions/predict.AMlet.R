predict.AMlet=function(X.old,fit_AMlet, X.new){
  ###################################################################
  ####library(wavethresh)
  ####source('softwavethresh.R')
  ####source('AMlet.R')
  ###################################################################
  P=ncol(X.old)
  if(ncol(X.new)!=P){
    print('PROBLEM:Not same nbr of predictors between training and test sets')
    out=NULL
  }
  
  mu.new.hat=matrix(NA,nrow = nrow(X.new),ncol = ncol(X.new))
  for (j in 1:P) {
    mu.new.hat[,j]=approx(X.old[,j],fit_AMlet$mu.hat[,j],xout = X.new[,j],rule=2)$y##careful!!
  }
  
  y.new=rowSums(mu.new.hat)+fit_AMlet$c.hat
  
  out=NULL
  out$mu.predict=mu.new.hat
  out$y.predict=y.new
  return(out)
}