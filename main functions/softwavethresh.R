softwavethresh<-function(y,X,lambda=NA,lambdatype='QUT',sigma=NA,p0=1,pen.father=F,
               filter.number=1,family="DaubExPhase", bc="periodic",
               alpha=0.05,M_MC=1000,M_GEV=100,
               type='GEV',com.sure=F,com.sqrtrss=F){
  N=length(y)
  J=log2(p0)
  
  Xorder=apply(X,2,order)
  Xrank=apply(X,2,rank)
  y=y[Xorder]
  
  if(is.na(sigma)){
    ytHaar=wd(y,1,family="DaubExPhase")  
    sigma=mad(accessD(ytHaar,log(N,2)-1),center=0)
  }
  
  if(is.na(lambda)){
      if(lambdatype=='QUT'){lambda=lambdaQUT.softwavethresh(X=X,p0=p0,pen.father=pen.father,sigma=sigma,filter.number=filter.number,family=family,bc= bc,
                                                            alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)}
      if(lambdatype=='Universal'){lambda=sigma*sqrt(2*log(N))}
      
  }
  
  yt=wd(y,filter.number=filter.number, family=family, bc=bc)
  yt.vec=rep(NA,N)
  yt.vec[1:p0]=accessC(yt,J)
  ii=p0
  for(level in J:(yt$nlevels-1)){
    lj=2^level
    yt.vec[(ii+1):(ii+lj)]=accessD(yt,level)
    ii=ii+lj
  }
  
  alpha.hat=rep(0,N)
  
  pos=which(yt.vec>lambda)
  neg=which(yt.vec<(-lambda))
  alpha.hat[pos]=yt.vec[pos]-lambda
  alpha.hat[neg]=yt.vec[neg]+lambda
  
  
  if(pen.father==F){
    alpha.hat[1:p0]=yt.vec[1:p0]
  }
  
  ytthresh=yt
  ytthresh=putC(ytthresh,J,alpha.hat[1:p0])
  ii=p0
  for(level in J:(ytthresh$nlevels-1)){
    lj=2^level
    ytthresh=putD(ytthresh, level,alpha.hat[(ii+1):(ii+lj)])
    ii=ii+lj
  }
  
  
  mu.hat=wr(ytthresh)[Xrank]
  
  out=NULL
  out$sigma=sigma
  out$lambdaQUT=lambda
  out$alpha.hat=alpha.hat
  out$mu.hat=mu.hat
  
  if(com.sure==T){
    sure.term.2=sum((yt.vec-alpha.hat)^2)
    out$sure.term.2=sure.term.2
    if(pen.father==F){sure.term.3.trace=p0+length(pos)+length(neg)}
    if(pen.father==T){sure.term.3.trace=length(pos)+length(neg)}
    out$sure.term.3.trace=sure.term.3.trace
    out$SURE=-N*sigma^2+sure.term.2+2*sigma^2*out$sure.term.3.trace
  }
  
  if(com.sqrtrss==T){
   out$sqrtrss=sqrt(mean((yt.vec-alpha.hat)^2))
  }
  
  return(out)
}
