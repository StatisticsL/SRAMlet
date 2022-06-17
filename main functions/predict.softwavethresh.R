predict.softwavethresh=function(X.old,fit_softwavethresh, X.new){
  ###################################################################
  ###################################################################
  if(ncol(X.new)!=1){
    print('PROBLEM:dimension has to be 1!')
    y.new=NULL
  }else{
    mu.hat=fit_softwavethresh$mu.hat
    y.new=approx(X.old,mu.hat,xout = X.new,rule=2)$y##careful!!
  }
  return(y.new)
}