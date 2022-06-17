lambdaQUT.softwavethresh<-function(X,p0=1,pen.father=F,sigma,filter.number=1,family="DaubExPhase", bc="periodic",
                                   alpha=0.05,M_MC=1000,M_GEV=100,type='GEV'){
  ##################################################################
  ###if pen.father=T, p0=0; if pen.father=F, p0=p0
  ##################################################################
  N=dim(X)[1]
  Xorder=apply(X,2,order)
  J=log2(p0)######check!!!!!
  
  numfct=function(z,xorder,pen.father,filter.number,family, bc){
    n=length(z)
    yt=wd(z[xorder],filter.number=filter.number,family=family, bc=bc)
    yt.vec=rep(NA,n)
    yt.vec[1:p0]=accessC(yt,J)
    yt.vec[(p0+1):n]=yt$D[1:(n-p0)]
    if(pen.father==F){return(max(abs(yt.vec[(p0+1):n])))}
    if(pen.father==T){return(max(abs(yt.vec)))}
  }

  if(type=='MC'){
    z=matrix(rnorm(N*M_MC,mean=0,sd=sigma), nrow=N, ncol=M_MC)
    z=t(t(z)-apply(z,2,mean))##??
    lambdas=apply(z,2,numfct,xorder=Xorder,pen.father=pen.father,filter.number=filter.number,family=family, bc=bc)
    lambdaQUT=quantile(lambdas,1-alpha)
  }
  if(type=='GEV'){
    z=matrix(rnorm(N*M_GEV,mean=0,sd=sigma), nrow=N, ncol=M_GEV)
    z=t(t(z)-apply(z,2,mean))##??
    lambdas=apply(z,2,numfct,xorder=Xorder,pen.father=pen.father,filter.number=filter.number,family=family, bc=bc)
    
    GEVfit <- gev.fit(lambdas,show=FALSE)
    GEVpars <- GEVfit$mle
    lambdaQUT <- Q(1-alpha, mu = GEVpars[1], sigma = GEVpars[2], xi = GEVpars[3])
  }
  return(lambdaQUT)
}