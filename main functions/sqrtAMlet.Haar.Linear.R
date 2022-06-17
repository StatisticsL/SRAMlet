sqrtAMlet.Haar.Linear<-function(y,X,lambda=NA,p0=1,filter.number=1,family="DaubExPhase",bc="periodic",
                           alpha=0.05,M_MC=1000,M_GEV=100,type='GEV',com.sqrtrss=F,
                           conv.thresh=1.e-10,max.iter=ncol(X)*10){
  N=length(y)
  P=ncol(X)
  com.sqrtrss=F
  
  if(P==1){print('Please use sqrtwavethresh function!');out=NULL;}
  if(P>1){##
    
    l2norm=sqrt(apply(X^2,2,sum))
    Xl2=X/matrix(rep(l2norm,N),N,P,byrow = TRUE)
    
    if(is.na(lambda)){
      lambdaSeris=lambdaQUT.sqrtAMlet.Haar.Linear(X=X,filter.number=filter.number,family=family, bc=bc,
                                                        alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
      lambda=lambdaSeris$lambdaQUT
    }
    lambda.SH=lambdaSeris$lambda.QUT.HS
    lambda.L=lambdaSeris$lambda.QUT.Linear
    
    cc=mean(y)
    muhat=matrix(0, N, 3*P)
    beta.hat=matrix(0,1,P)##linear part
    alpha.hat=matrix(0,N,P)###smooth part
    gamma.hat=matrix(0,N,P)####haar part
    res=y-rowSums(muhat)-cc
    
    First=T
    continue=T
    indexes=seq(1,3*P)
    index.deleted=3*P+1
    iter=0
    
    while (continue&(iter<max.iter)) {
     
      muhat.orig=muhat
      
      for (p_index in indexes[-index.deleted]) {
        resp=res+muhat[,p_index]
        if(p_index <= P){
          outp=sqrtLinear(resp*lambda.L/lambda.SH,as.matrix(Xl2[,p_index]),lambda*lambda.L/lambda.SH,com.sqrtrss=com.sqrtrss)
          beta.hat[1,p_index]=outp$beta.hat
          muhat[,p_index]=outp$mu.hat*lambda.SH/lambda.L
          # muhat[,p_index]=outp$mu.hat
        }else if(p_index <= 2*P){
          outp=sqrtwavethresh(y=resp,as.matrix(X[,(p_index-P)]),lambda=lambda,sigma=NA,p0=p0,pen.father=T,
                              filter.number=filter.number, family= family,bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
          alpha.hat[,(p_index-P)]=outp$alpha.hat
          muhat[,p_index]=outp$mu.hat
        }else{
          outp=sqrtwavethresh(y=resp,as.matrix(X[,(p_index-2*P)]),lambda=lambda,sigma=NA,p0=p0,pen.father=T,
                              filter.number=1, family= "DaubExPhase",bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
          gamma.hat[,(p_index-2*P)]=outp$alpha.hat
          muhat[,p_index]=outp$mu.hat
        }
        # if(com.sqrtrss==T){sqrtrss.vec[p_index]=outp$sqrtrss}
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
        index.deleted=3*P+1
      }
    }
    
    #################################################
    #################################################
    ###################?????????????????????
    # for (p_index in 1:P) {
    #   beta.hat[1,p_index]=beta.hat[1,p_index]/l2norm[p_index]
    #   }
    ###########------results-----------------------------------------
    out=NULL
    out$c.hat=cc
    out$alpha.hat=alpha.hat
    out$beta.hat=beta.hat
    out$gamma.hat=gamma.hat
    muhat.aa=muhat[,1:P];muhat.bb=muhat[,(P+1):(2*P)];muhat.cc=muhat[,(2*P+1):(3*P)];
    out$mu.hat=muhat.aa+muhat.bb+muhat.cc
    out$lambdaQUT=lambda
    out$index.selected=which((colSums(abs(alpha.hat)+abs(gamma.hat))+abs(beta.hat))!=0)
  }
  return(out)
}

