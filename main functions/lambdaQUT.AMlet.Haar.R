lambdaQUT.AMlet.Haar<-function(X,sigma,filter.number=1,family="DaubExPhase", bc="periodic",
                          alpha=0.05,M_MC=1000,M_GEV=100,type='GEV'){
  ##################################################################
  ###p0=1
  ##################################################################
  N=dim(X)[1]
  Xorder=apply(X,2,order)
  
  numj=function(z,xorder,filter.number,family, bc){
    n=length(z)
    yt=wd(z[xorder],filter.number=filter.number,family=family, bc=bc)
    yt.vec=rep(NA,n)
    yt.vec[1]=accessC(yt,0)
    yt.vec[2:n]=yt$D[1:(n-1)]
    return(max(abs(yt.vec)))
  }
  
  
  numfct=function(z,Xorder,filter.number,family, bc){
    n=length(z)
    num.vec=apply(Xorder,2,numj,z=z,filter.number=filter.number,family=family,bc=bc)
    return(max(num.vec))
  }
  
  if(type=='MC'){
    z=matrix(rnorm(N*M_MC,mean = 0,sd=sigma), nrow=N, ncol=M_MC)
    z=t(t(z)-apply(z,2,mean))##??
    lambdas1=apply(z,2,numfct,Xorder=Xorder,filter.number=filter.number,family=family, bc=bc)
    lambdas2=apply(z,2,numfct,Xorder=Xorder,filter.number=1,family="DaubExPhase", bc=bc)
    lambdas=apply(matrix(c(lambdas1,lambdas2),2,M_MC,byrow = TRUE),2,max)
    lambdaQUT=quantile(lambdas,1-alpha)
  }
  if(type=='GEV'){
    z=matrix(rnorm(N*M_GEV,mean = 0,sd=sigma), nrow=N, ncol=M_GEV)
    z=t(t(z)-apply(z,2,mean))##??
    lambdas1=apply(z,2,numfct,Xorder=Xorder,filter.number=filter.number,family=family, bc=bc)
    lambdas2=apply(z,2,numfct,Xorder=Xorder,filter.number=1,family="DaubExPhase", bc=bc)
    lambdas=apply(matrix(c(lambdas1,lambdas2),2,M_GEV,byrow = TRUE),2,max)
    
    GEVfit <- gev.fit(lambdas,show=FALSE)
    GEVpars <- GEVfit$mle
    lambdaQUT <- Q(1-alpha, mu = GEVpars[1], sigma = GEVpars[2], xi = GEVpars[3])
  }
  return(lambdaQUT)
}