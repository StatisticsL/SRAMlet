predict.sqrtwavethresh=function(X.old,fit_sqrtwavethresh, X.new){
  ###################################################################
  ###################################################################
  if(ncol(X.new)!=1){
    print('PROBLEM:dimension has to be 1!')
    y.new=NULL
  }else{
    mu.hat=fit_sqrtwavethresh$mu.hat
    y.new=approx(X.old,mu.hat,xout = X.new,rule=2)$y##careful!!
  }
  return(y.new)
}