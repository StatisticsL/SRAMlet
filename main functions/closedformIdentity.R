closedformIdentity=function(z,lambda,p0,com.sure=F,com.sqrtrss=F){
  N=length(z)
  if(p0==0){dataXi=z}
  if(p0>0){
    dataGamma=z[1:p0]
    dataXi=z[(p0+1):N]
  }
  
  lambdaXi=1/sqrt(N-p0-sum(dataXi==0))
  lambdaZero=max(abs(dataXi))/sqrt(sum(dataXi^2))
  
  if(lambda<=lambdaXi){
    w=dataXi
    if(com.sure==T){sure.term.3.trace=N}
  }else if(lambda >= lambdaZero){
      w=rep(0,N-p0)
      if(com.sure==T){sure.term.3.trace=p0}
  }else{
      sdataXi=sort(dataXi^2)###increasing
      pu=cumsum(sdataXi)[-(N-p0)]/(1-lambda^2*seq(N-p0-1,1,-1))
      i0=sum(pu<=0)+sum(pu==Inf)+1
      ind=i0:(N-p0-1)
      varphis=lambda*sqrt(pu[ind])
      j0=which((n0(dataXi,varphis)==ind)==T)
      if(length(j0)==1){
        varphis0=varphis[j0]
        w=rep(0,(N-p0))
        wpos=(dataXi>varphis0)
        wneg=(dataXi< -varphis0)
        w[wpos]=dataXi[wpos]-varphis0
        w[wneg]=dataXi[wneg]+varphis0
        if(com.sure==T){
          wnonzero=which(w!=0)
          lwn=length(wnonzero)
          sure.term.3.trace=p0+lwn
          
          # Jxi=matrix(0,nrow = (N-p0),ncol = (N-p0))
          # Jxi[wnonzero,wnonzero]=-lambda^2
          # if(lwn!=1){
          #   diag(Jxi[wnonzero,wnonzero])=1-lambda^2
          # }else{
          #   Jxi[wnonzero,wnonzero]=1-lambda^2
          # }
          # 
          # Jwinv=diag(N-p0)
          # Jwinv[wnonzero,wnonzero]=solve(Jxi[wnonzero,wnonzero])
          # Jwinv[wnonzero,wnonzero]=Jwinv[wnonzero,wnonzero]-lambda^2
          # Jwinv[wnonzero,wnonzero]=solve(Jwinv[wnonzero,wnonzero])
          # sure.term.3.trace=p0+sum(diag(Jwinv%*%Jxi))
        }
      }else if(length(j0)<1){
        print('Problem: can find j!')
        w=rep(NA,N-p0)
      }else{
        print('Problem: j is not unique!')
        w=rep(NA,N-p0)
      }
  }
  if(p0==0){mu=w}
  if(p0>0){mu=c(dataGamma,w)}
  
  if(com.sqrtrss==T){sqrtrss=sqrt(mean((z-mu)^2))}
  
  out=NULL
  out$mu=mu
  if(com.sure==T){
    out$sure.term.3.trace=sure.term.3.trace
    out$sure.term.2=sum((z-mu)^2)
    }
  if(com.sqrtrss==T){out$sqrtrss=sqrtrss}
  
  return(out)
}
