sqrtLinear<-function(y,x,lambda,com.sqrtrss=F){
  ###x l2-normalised
  lambda0=abs(sum(x*y))/sqrt(sum(y^2))
  if(lambda >= lambda0){
    betavalue=0
  }else{
    diffe=sum(y^2)-(sum(x*y))^2
    if(abs(diffe)<=1e-8){ diffe = 0}
    varphi= lambda/sqrt(1-lambda^2)*sqrt(diffe)
    if(sum(x*y)>varphi){betavalue=sum(x*y)-varphi}
    if(sum(x*y)< -varphi){betavalue=sum(x*y)+varphi}
  }
  
  out=NULL
  out$beta.hat=betavalue
  out$mu.hat=betavalue*x
  if(com.sqrtrss==TRUE){
    out$sqrtrss=sqrt(mean((y-betavalue*x)^2))
    }
  return(out)
}