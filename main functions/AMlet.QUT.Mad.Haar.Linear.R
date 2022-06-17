AMlet.QUT.Mad.Haar.Linear<-function(y,X,lambda=NA,sigma=NA,p0=1,filter.number=1,family="DaubExPhase",bc="periodic",
                               alpha=0.05,M_MC=1000,M_GEV=100,type='GEV',com.sqrtrss=F,
                               conv.thresh=1.e-10, max.iter=ncol(X)*10){
  N=length(y)
  P=ncol(X)
  com.sqrtrss=F
  
  l2norm=sqrt(apply(X^2,2,sum))
  Xl2=X/matrix(rep(l2norm,N),N,P,byrow = TRUE)
  
  if(P==1){print('Please use softwavethresh function!');out=NULL;}
  if(P>1){
    if(is.na(lambda)&is.na(sigma)){
      lambdaSeris=lambdaQUT.AMlet.Haar.Linear(X=X,sigma=1,filter.number=filter.number,family=family, bc=bc,
                                     alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
      lambda0=lambdaSeris$lambdaQUT
      lambda0.HS=lambdaSeris$lambda.QUT.HS
      lambda0.L=lambdaSeris$lambda.QUT.Linear
      
      cc=mean(y)
      muhat=matrix(0, N,3*P)###linear,smooth,haar
      res=y-rowSums(muhat)-cc
      sigma.hat=rep(0,3*P)
      
      for (p_index in 1:(3*P)) {
        resp=res+muhat[,p_index]
        if(p_index <= P){##haar
          outp=softwavethresh(y=resp,as.matrix(X[,p_index]),lambda=NA,sigma=NA,p0=p0,pen.father=T,
                              filter.number=1,family="DaubExPhase",bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
          sigma.hat[p_index]=outp$sigma
        }else if((p_index > P)&(p_index <= (2*P))){##smooth
          outp=softwavethresh(y=resp,as.matrix(X[,(p_index-P)]),lambda=NA,sigma=NA,p0=p0,pen.father=T,
                              filter.number=filter.number,family=family,bc=bc,
                              alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
          sigma.hat[p_index]=outp$sigma
        }else{##
          ytHaar=wd(resp,1,family="DaubExPhase")  
          sigma.hat[p_index]=mad(accessD(ytHaar,log(N,2)-1),center=0)
        }
        
        res=resp-muhat[,p_index]
      }
      
      cc=mean(y)
      muhat=matrix(0, N, 3*P)
      alpha.hat=matrix(0,N,P)##smooth 
      beta.hat=matrix(0,1,P)##linear
      gamma.hat=matrix(0,N,P)##haar
      res=y-rowSums(muhat)-cc
      
      First=T
      continue=T
      indexes=seq(1,3*P)
      index.deleted=3*P+1
      iter=0
      
      while (continue&(iter<max.iter)) {
        
        sigma=sum(sigma.hat)/(3*P)
        lambda=lambda0*sigma
        lambda.HS=lambda0.HS*sigma
        lambda.L=lambda0.L*sigma
        
        muhat.orig=muhat
        
        for (p_index in indexes[-index.deleted]) {
          resp=res+muhat[,p_index]
          if(p_index<= P){
            outp=softwavethresh(y=resp,as.matrix(X[,p_index]),lambda=lambda,sigma=sigma,p0=p0,pen.father=T,
                                filter.number=1,family="DaubExPhase",bc=bc,
                                alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
            gamma.hat[,p_index]=outp$alpha.hat
            muhat[,p_index]=outp$mu.hat
            sigma.hat[p_index]=outp$sigma
          }else if((p_index > P)&(p_index <= (2*P))){
            outp=softwavethresh(y=resp,as.matrix(X[,(p_index-P)]),lambda=NA,sigma=NA,p0=p0,pen.father=T,
                                filter.number=filter.number,family=family,bc=bc,
                                alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
            alpha.hat[,(p_index-P)]=outp$alpha.hat
            sigma.hat[p_index]=outp$sigma
            muhat[,p_index]=outp$mu.hat
          }else{
            beta.hat[1,(p_index-2*P)]=softfunction(sum(resp*lambda.L/lambda.HS*as.matrix( Xl2[,(p_index-2*P)])),lambda*lambda.L/lambda.HS)/sum(as.matrix( Xl2[,(p_index-2*P)])^2)
            muhat[,p_index]=beta.hat[1,(p_index-2*P)]*as.matrix( Xl2[,(p_index-2*P)])*lambda.HS/lambda.L
            # muhat[,p_index]=beta.hat[1,(p_index-2*P)]*as.matrix( Xl2[,(p_index-2*P)])
            ytHaar=wd(resp,1,family="DaubExPhase")  
            sigma.hat[p_index]=mad(accessD(ytHaar,log(N,2)-1),center=0)
          }
          
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
    }else{
      if(is.na(lambda)){
        lambdaSeris=lambdaQUT.AMlet.Haar.Linear(X=X,sigma=sigma,filter.number=filter.number,family=family, bc=bc,
                                                      alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type)
        lambda=lambdaSeris$lambdaQUT
        lambda.HS=lambdaSeris$lambda.QUT.HS
        lambda.L=lambdaSeris$lambda.QUT.Linear
        }
      cc=mean(y)
      muhat=matrix(0, N, 3*P)
      alpha.hat=matrix(0,N,P)
      beta.hat=matrix(0,1,P)
      gamma.hat=matrix(0,N,P)
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
          if(p_index<= P){
            outp=softwavethresh(y=resp,as.matrix(X[,p_index]),lambda=lambda,sigma=sigma,p0=p0,pen.father=T,
                                filter.number=1,family="DaubExPhase",bc=bc,
                                alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
            gamma.hat[,p_index]=outp$alpha.hat
            muhat[,p_index]=outp$mu.hat
          }else if((p_index > P)&(p_index <= (2*P))){
            outp=softwavethresh(y=resp,as.matrix(X[,(p_index-P)]),lambda=NA,sigma=NA,p0=p0,pen.father=T,
                                filter.number=filter.number,family=family,bc=bc,
                                alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,type=type,com.sure=F,com.sqrtrss=com.sqrtrss)
            alpha.hat[,(p_index-P)]=outp$alpha.hat
            muhat[,p_index]=outp$mu.hat
          }else{
            beta.hat[1,(p_index-2*P)]=softfunction(sum(resp*lambda.L/lambda.HS * as.matrix( Xl2[,(p_index-2*P)])),lambda*lambda.L/lambda.HS)/sum(as.matrix( Xl2[,(p_index-2*P)])^2)
            muhat[,p_index]=beta.hat[1,(p_index-2*P)]*as.matrix( Xl2[,(p_index-2*P)])*lambda.HS/lambda.L
            # muhat[,p_index]=beta.hat[1,(p_index-2*P)]*as.matrix( Xl2[,(p_index-2*P)])
          }
          
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
    }
    
    #################################################
    #################################################
    #################################################
    # for (p_index in 1:P) {
    #   beta.hat[1,p_index]=beta.hat[1,p_index]/l2norm[p_index]
    # }
    
    ###########------results-----------------------------------------
    out=NULL
    out$c.hat=cc
    out$alpha.hat=alpha.hat
    out$beta.hat=beta.hat
    out$gamma.hat=gamma.hat
    muhat.aa=muhat[,1:P];muhat.bb=muhat[,(P+1):(2*P)];muhat.cc=muhat[,(2*P+1):(3*P)];
    out$mu.hat=muhat.aa+muhat.bb+muhat.cc
    out$lambdaQUT=lambda
    out$sigma=sigma
    out$index.selected=which((colSums(abs(alpha.hat)+abs(gamma.hat))+abs(beta.hat))!=0)
  }
  return(out)
}

