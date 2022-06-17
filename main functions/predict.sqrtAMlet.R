predict.sqrtAMlet=function(X.old,fit_sqrtAMlet, X.new){
  ###################################################################
  ####library(wavethresh)
  ####source('sqrtwavethresh.R')
  ####source('sqrtAMlet.R')
  ###################################################################
  P=ncol(X.old)
  if(ncol(X.new)!=P){
    print('PROBLEM:Not same nbr of predictors between training and test sets')
    out=NULL
  }else{
    mu.new.hat=matrix(NA,nrow = nrow(X.new),ncol = ncol(X.new))
    mu.hat=fit_sqrtAMlet$mu.hat
    for (j in 1:P) {
      mu.new.hat[,j]=approx(X.old[,j],mu.hat[,j],xout = X.new[,j],rule=2)$y##careful!!
    }
    
    y.new=rowSums(mu.new.hat)+fit_sqrtAMlet$c.hat
    out=NULL
    out$mu.predict=mu.new.hat
    out$y.predict=y.new
  }
  return(out)
}