SRAMlet<-function(y,X,lambda=NA,p0=1,addHaar=T,addLinear=T,
                    filter.number=1,family="DaubExPhase",
                    alpha=0.05,M_MC=1000,M_GEV=100,type='GEV',
                    conv.thresh=1.e-10, max.iter=ncol(X)*10,com.sqrtrss=F){
  ###com.sqrtrss is just for only smooth part, not for adding other part
  bc="periodic";
  if((addHaar==T) & (addLinear==T)){
    out=sqrtAMlet.Haar.Linear(y=y,X=X,lambda=lambda,p0=p0,filter.number=filter.number,family=family,bc=bc,
                                        alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sqrtrss=com.sqrtrss,
                                        conv.thresh=conv.thresh,max.iter=max.iter)
  }else if((addHaar==T) &(addLinear==F)){
    out=sqrtAMlet.Haar(y=y,X=X,lambda=lambda,p0=p0,filter.number=filter.number,family=family,bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sqrtrss=com.sqrtrss,
                              conv.thresh=conv.thresh,max.iter=max.iter)
  }else if((addHaar==F) & (addLinear==T)){
    out=sqrtAMlet.Linear(y=y,X=X,lambda=lambda,p0=p0,filter.number=filter.number,family=family,bc=bc,
                       alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sqrtrss=com.sqrtrss,
                       conv.thresh=conv.thresh,max.iter=max.iter)
  }else{
    out=sqrtAMlet(y=y,X=X,lambda=lambda,p0=p0,filter.number=filter.number,family=family,bc=bc,
                         alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sqrtrss=com.sqrtrss,
                         conv.thresh=conv.thresh,max.iter=max.iter)
  }
  
  return(out)
}