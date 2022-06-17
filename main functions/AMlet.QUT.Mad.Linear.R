AMlet.QUT.Mad.Linear<-function(y,X,lambda=NA,sigma=NA,p0=1,filter.number=1,family="DaubExPhase",bc="periodic",
                             alpha=0.05,M_MC=1000,M_GEV=100,type='GEV',com.sqrtrss=F,
                             conv.thresh=1.e-10, max.iter=ncol(X)*10){
  N=length(y)
  P=ncol(X)
  
  l2norm=sqrt(apply(X^2,2,sum))
  Xl2=X/matrix(rep(l2norm,N),N,P,byrow = TRUE)
  
  if(P==1){print('Please use softwavethresh function!');out=NULL;}
  if(P>1){
    if(is.na(lambda)&is.na(sigma)){
      lambdaSeris=lambdaQUT.AMlet.Linear(X=X,sigma=1,filter.number=filter.number,family=family, bc=bc,
                                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
      lambda0=lambdaSeris$lambdaQUT
      lambda0.S=lambdaSeris$lambda.QUT.S
      lambda0.L=lambdaSeris$lambda.QUT.Linear
      
      cc=mean(y)
      muhat=matrix(0, N,2*P)
      res=y-rowSums(muhat)-cc
      sigma.hat=rep(0,2*P)
      
      for (p_index in 1:(2*P)) {
        resp=res+muhat[,p_index]
        if(p_index <= P){
          outp=softwavethresh(y=resp,as.matrix(X[,p_index]),lambda=NA,sigma=NA,p0=p0,pen.father=T,
                              filter.number=filter.number,family=family,bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
          sigma.hat[p_index]=outp$sigma
        }else{
          ytHaar=wd(resp,1,family="DaubExPhase")  
          sigma.hat[p_index]=mad(accessD(ytHaar,log(N,2)-1),center=0)
        }
        res=resp-muhat[,p_index]
      }
      
      cc=mean(y)
      muhat=matrix(0, N, 2*P)
      alpha.hat=matrix(0,N,P)
      beta.hat=matrix(0,1,P)
      res=y-rowSums(muhat)-cc
      
      First=T
      continue=T
      indexes=seq(1,2*P)
      index.deleted=2*P+1
      iter=0
      if(com.sqrtrss==T){sqrtrss.vec=rep(NA,2*P)}
      
      while (continue&(iter<max.iter)) {
        
        sigma=sum(sigma.hat)/(2*P)
        lambda=lambda0*sigma
        lambda.S=lambda0.S*sigma
        lambda.L=lambda0.L*sigma
        
        muhat.orig=muhat
        
        for (p_index in indexes[-index.deleted]) {
          resp=res+muhat[,p_index]
          if(p_index<= P){
            outp=softwavethresh(y=resp,as.matrix(X[,p_index]),lambda=lambda,sigma=sigma,p0=p0,pen.father=T,
                                filter.number=filter.number,family=family,bc=bc,
                                alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
            alpha.hat[,p_index]=outp$alpha.hat
            muhat[,p_index]=outp$mu.hat
            sigma.hat[p_index]=outp$sigma
          }else{
            beta.hat[1,(p_index-P)]=softfunction(sum(resp*lambda.L/lambda.S*as.matrix(Xl2[,(p_index-P)])),lambda*lambda.L/lambda.S)/sum(as.matrix(Xl2[,(p_index-P)])^2)
            muhat[,p_index]=beta.hat[1,(p_index-P)]*as.matrix(Xl2[,(p_index-P)])*lambda.S/lambda.L
            ytHaar=wd(resp,1,family="DaubExPhase")  
            sigma.hat[p_index]=mad(accessD(ytHaar,log(N,2)-1),center=0)
          }
          
          if(com.sqrtrss==T){sqrtrss.vec[p_index]=outp$sqrtrss}
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
    }else{
      if(is.na(lambda)){
        lambdaSeris=lambdaQUT.AMlet.Linear(X=X,sigma=sigma,filter.number=filter.number,family=family, bc=bc,
                                           alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
        lambda=lambdaSeris$lambdaQUT
        lambda.S=lambdaSeris$lambda.QUT.S
        lambda.L=lambdaSeris$lambda.QUT.Linear}
      cc=mean(y)
      muhat=matrix(0, N, 2*P)
      alpha.hat=matrix(0,N,P)
      beta.hat=matrix(0,1,P)
      res=y-rowSums(muhat)-cc
      
      First=T
      continue=T
      indexes=seq(1,2*P)
      index.deleted=2*P+1
      iter=0
      
      if(com.sqrtrss==T){sqrtrss.vec=rep(NA,2*P)}
      
      while (continue&(iter<max.iter)) {
        
        muhat.orig=muhat
        
        for (p_index in indexes[-index.deleted]) {
          resp=res+muhat[,p_index]
          if(p_index <= P){
            outp=softwavethresh(y=resp,as.matrix(X[,p_index]),lambda=lambda,sigma=sigma,p0=p0,pen.father=T,
                                filter.number=filter.number,family=family,bc=bc,
                                alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
            alpha.hat[,p_index]=outp$alpha.hat
            muhat[,p_index]=outp$mu.hat
          }else{
            beta.hat[1,(p_index-P)]=softfunction(sum(resp*lambda.L/lambda.S*as.matrix(Xl2[,(p_index-P)])),lambdalambda.L/lambda.S)/sum(as.matrix(Xl2[,(p_index-P)])^2)
            muhat[,p_index]=beta.hat[1,(p_index-P)]*as.matrix(Xl2[,(p_index-P)])*lambda.S/lambda.L
          }
          
          if(com.sqrtrss==T){sqrtrss.vec[p_index]=outp$sqrtrss}
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
    }
    
    #################################################
    #################################################
    #################################################
    # for (p_index in 1:P) {
    #   beta.hat[1,p_index]=beta.hat[1,p_index]/l2norm[p_index]
    # 
    # }
    ###########------results-----------------------------------------
    out=NULL
    out$c.hat=cc
    out$alpha.hat=alpha.hat
    out$beta.hat=beta.hat
    muhat.aa=muhat[,1:P];muhat.bb=muhat[,(P+1):(2*P)]
    out$mu.hat=muhat.aa+muhat.bb
    out$lambdaQUT=lambda
    out$sigma=sigma
    out$index.selected=which((colSums(abs(alpha.hat))+abs(beta.hat))!=0)
  }
  return(out)
}

