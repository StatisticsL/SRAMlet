lambdaQUT.sqrtAMlet.Haar.Linear<-function(X,filter.number=1,family="DaubExPhase", bc="periodic",
                                     alpha=0.05,M_MC=1000,M_GEV=100,type='GEV'){
  ##################################################################
  ###p0=1
  ##################################################################
  N=dim(X)[1]
  P=ncol(X)
  l2norm=sqrt(apply(X^2,2,sum))
  Xl2=X/matrix(rep(l2norm,N),N,P,byrow = TRUE)
  
  lambda.QUT.HS=lambdaQUT.sqrtAMlet.Haar(X=X,filter.number=filter.number,family=family, bc=bc,
                                                alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
  
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
    z=matrix(rnorm(N*M_MC), nrow=N, ncol=M_MC)
    z=t(t(z)-apply(z,2,mean))##??
    den=sqrt(apply(z^2,2,sum))
    num=apply(z,2,numfct,Xorder=Xorder,filter.number=filter.number,family=family, bc=bc)
    lambdas1=num/den
    num2=t(Xl2)%*%z
    lambdas2=apply(abs(num2),2,max)/den
    lambdas.QUT.linear=quantile(lambdas2,1-alpha)
    num2=t(Xl2*lambda.QUT.HS/lambdas.QUT.linear)%*%z
    lambdas2=apply(abs(num2),2,max)/den
    num3=apply(z,2,numfct,Xorder=Xorder,filter.number=1,family="DaubExPhase", bc=bc)
    lambdas3=num3/den
    lambdas=apply(matrix(c(lambdas1,lambdas2,lambdas3),3,M_MC,byrow = TRUE),2,max)
    lambdaQUT=quantile(lambdas,1-alpha)
  }
  if(type=='GEV'){
    z=matrix(rnorm(N*M_GEV), nrow=N, ncol=M_GEV)
    z=t(t(z)-apply(z,2,mean))##??
    # Smooth
    den=sqrt(apply(z^2,2,sum))
    num=apply(z,2,numfct,Xorder=Xorder,filter.number=filter.number,family=family, bc=bc)
    lambdas1=num/den
    # Linear
    num2=t(Xl2)%*%z
    lambdas2=apply(abs(num2),2,max)/den
    GEVfit.linear <- gev.fit(lambdas2,show=FALSE)
    GEVpars.linear<- GEVfit.linear$mle
    lambdas.QUT.linear<- Q(1-alpha, mu = GEVpars.linear[1], sigma = GEVpars.linear[2], xi = GEVpars.linear[3])
    num2=t(Xl2*lambda.QUT.HS/lambdas.QUT.linear)%*%z
    lambdas2=apply(abs(num2),2,max)/den
    # Haar
    num3=apply(z,2,numfct,Xorder=Xorder,filter.number=1,family="DaubExPhase", bc=bc)
    lambdas3=num3/den
    lambdas=apply(matrix(c(lambdas1,lambdas2,lambdas3),3,M_MC,byrow = TRUE),2,max)
    
    GEVfit <- gev.fit(lambdas,show=FALSE)
    GEVpars <- GEVfit$mle
    lambdaQUT <- Q(1-alpha, mu = GEVpars[1], sigma = GEVpars[2], xi = GEVpars[3])
  }
  out=NULL;
  out$lambda.QUT.HS=lambda.QUT.HS
  out$lambda.QUT.Linear=lambdas.QUT.linear
  out$lambdaQUT=lambdaQUT
  return(out)
}