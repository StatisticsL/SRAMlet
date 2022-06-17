sqrtAMlet.Haar<-function(y,X,lambda=NA,p0=1,filter.number=1,family="DaubExPhase",bc="periodic",
                    alpha=0.05,M_MC=1000,M_GEV=100,type='GEV',com.sqrtrss=F,
                    conv.thresh=1.e-10, max.iter=ncol(X)*10){
  N=length(y)
  P=ncol(X)
  com.sqrtrss=F
  if(P==1){print('Please use sqrtwavethresh function!');out=NULL;}
  if(P>1){###problem: lambda need to be changed if adding Haar!!!!
    if(is.na(lambda)){lambda=lambdaQUT.sqrtAMlet.Haar(X,filter.number=filter.number,family=family, bc=bc,
                      alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)}
    cc=mean(y)
    muhat=matrix(0,N,2*P)
    alpha.hat=matrix(0,N,2*P)
    res=y-rowSums(muhat)-cc
    
    
    First=T
    continue=T
    indexes=seq(1,2*P)
    index.deleted=2*P+1
    iter=0
    
    # if(com.sqrtrss==T){sqrtrss.vec=rep(NA,2*P)}
    
    while (continue&(iter<max.iter)) {
      
      muhat.orig=muhat
      
      for (p_index in indexes[-index.deleted]) {
        resp=res+muhat[,p_index]
        if(p_index <= P){
          outp=sqrtwavethresh(y=resp,as.matrix(X[,p_index]),lambda=lambda,sigma=NA,p0=p0,pen.father=T,
                              filter.number=1, family="DaubExPhase",bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
        }else{
          outp=sqrtwavethresh(y=resp,as.matrix(X[,(p_index-P)]),lambda=lambda,sigma=NA,p0=p0,pen.father=T,
                              filter.number=filter.number, family= family,bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
        }
        # if(com.sqrtrss==T){sqrtrss.vec[p_index]=outp$sqrtrss}
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
        index.deleted=2*P+1
      }
    }
    
    ###########------results-----------------------------------------
    out=NULL
    out$c.hat=cc
    out$alpha.hat=alpha.hat
    muhat.aa=muhat[,1:P];muhat.bb=muhat[,(P+1):(2*P)]
    out$mu.hat=muhat.aa+muhat.bb
    out$lambdaQUT=lambda
    # if(com.sqrtrss==T){
    #   out$sqrtrss.vec.2P=sqrtrss.vec
    #   sqrtrss.vec.aa=sqrtrss.vec[1:P];sqrtrss.vec.bb=sqrtrss.vec[(P+1):(2*P)]
    #   out$sqrtrss.vec=sqrt((sqrtrss.vec.aa^2+sqrtrss.vec.bb^2)/2)
    #   out$sqrtrss=sqrt(mean((y-rowSums(muhat))^2))
    # }
    alpha.hat.aa=abs(alpha.hat[,1:P]); alpha.hat.bb=abs(alpha.hat[,(P+1):(2*P)])
    out$index.selected=which(colSums(alpha.hat.aa+ alpha.hat.bb)!=0)
  }
  return(out)
}

