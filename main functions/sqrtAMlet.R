sqrtAMlet<-function(y,X,lambda=NA,p0=1,filter.number=1,family="DaubExPhase",bc="periodic",
                     alpha=0.05,M_MC=1000,M_GEV=100,type='GEV',com.sqrtrss=F,
                     conv.thresh=1.e-10, max.iter=ncol(X)*10){
  N=length(y)
  P=ncol(X)
  
  if(P==1){print('Please use sqrtwavethresh function!');out=NULL;}
  if(P>1){
    if(is.na(lambda)){lambda=lambdaQUT.sqrtAMlet(X=X,filter.number=filter.number,family=family, bc=bc,
                                                 alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)}
    
    cc=mean(y)
    muhat=matrix(0, N, P)
    alpha.hat=matrix(0,N, P)
    res=y-rowSums(muhat)-cc
    
    First=T
    continue=T
    indexes=seq(1,P)
    index.deleted=P+1
    iter=0
    
    if(com.sqrtrss==T){sqrtrss.vec=rep(NA,P)}
    
    while (continue&(iter<max.iter)) {
      
      muhat.orig=muhat
      
      for (p_index in indexes[-index.deleted]) {
        resp=res+muhat[,p_index]
        outp=sqrtwavethresh(y=resp,as.matrix(X[,p_index]),lambda=lambda,sigma=NA,p0=p0,pen.father=T,
                            filter.number=filter.number,family=family,bc=bc,
                            alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
        if(com.sqrtrss==T){sqrtrss.vec[p_index]=outp$sqrtrss}
        muhat[,p_index]=outp$mu.hat
        alpha.hat[,p_index]=outp$alpha.hat
        res=resp-muhat[,p_index]
        if(sum(abs(muhat[,p_index]))==0){index.deleted=c(p_index,index.deleted)}
      }
      
      cc.new=mean(cc+res)
      res=res+cc-cc.new
      cc=cc.new
      iter=iter+1
      
      rela.error=mean(abs(muhat.orig-muhat))/(1e-8+mean(abs(muhat)))
      continue=(rela.error>conv.thresh)
      if(((continue==F)|(iter>max.iter)) & (First==T)){
        continue=T
        First=F
        iter=0
        index.deleted=P+1
      }
    }
    
    ###########------results-----------------------------------------
    out=NULL
    out$c.hat=cc
    out$alpha.hat=alpha.hat
    out$mu.hat=muhat
    out$lambdaQUT=lambda
    if(com.sqrtrss==T){
      out$sqrtrss.vec=sqrtrss.vec
      out$sqrtrss=sqrt(mean((y-rowSums(muhat))^2))
    }
   out$index.selected=which(colSums(abs(alpha.hat))!=0)
  }
  return(out)
}
 
