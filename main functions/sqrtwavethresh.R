sqrtwavethresh<-function(y,X,lambda=NA,sigma=NA,p0=1,pen.father=F,
               filter.number=1,family="DaubExPhase",bc="periodic",
               alpha=0.05,M_MC=1000,M_GEV=100,
               type='GEV',com.sure=F,com.sqrtrss=F){
  ###############################################################
  N=length(y)
  P=dim(X)[2]
  Xorder=apply(X,2,order)
  Xrank=apply(X,2,rank)
  y=y[Xorder]
  J=log2(p0)
  
  if(P!=1){
    print('the dimension has to be 1!')
    out=NULL
  }
  if(P==1){
    ########------compute lambdaQUT---------------------
    if(is.na(lambda)){
        lambda=lambdaQUT.sqrtwavethresh(X=X,p0=p0,pen.father=pen.father,filter.number=filter.number,family=family, bc=bc,
                                        alpha= alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
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
    
    if(pen.father==T){results=closedformIdentity(z=yt.vec,lambda=lambda,p0=0,com.sure=com.sure,com.sqrtrss=com.sqrtrss)}
    if(pen.father==F){results=closedformIdentity(z=yt.vec,lambda=lambda,p0=p0,com.sure=com.sure,com.sqrtrss=com.sqrtrss)}
    
    alpha.hat=results$mu
    
    ytthresh=yt
    ytthresh=putC(ytthresh,J,alpha.hat[1:p0])
    ii=p0
    for(level in J:(yt$nlevels-1)){
      lj=2^level
      ytthresh=putD(ytthresh, level,alpha.hat[(ii+1):(ii+lj)])
      ii=ii+lj
    }
    
    mu.hat=wr(ytthresh)[Xrank]
    
    out=NULL
    out$alpha.hat=alpha.hat
    out$mu.hat=mu.hat
    out$lambdaQUT=lambda
    if(com.sure==T){
      if(is.na(sigma)){
        ytHaar=wd(y,1,family="DaubExPhase")  
        sigma=mad(accessD(ytHaar,log(N,2)-1),center=0)
      }
      out$sigma=sigma
      out$sure.term.2=results$sure.term.2
      out$sure.term.3.trace=results$sure.term.3.trace
      out$SURE=-N*sigma^2+results$sure.term.2+2*sigma^2*results$sure.term.3.trace
    }
    if(com.sqrtrss){
      out$sqrtrss=results$sqrtrss
    }
  }
  return(out)
}