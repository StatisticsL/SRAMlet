linearsolution<-function(y,X,com.sqrtrss=F){
  ###############################################################
    alpha.hat=y/X
    mu.hat=alpha.hat*X
    
    out=NULL
    out$alpha.hat=alpha.hat
    out$mu.hat=mu.hat
    if(com.sqrtrss){
      out$sqrtrss=sqrt(mean((y-mu.hat)^2))
    }
  return(out)
}