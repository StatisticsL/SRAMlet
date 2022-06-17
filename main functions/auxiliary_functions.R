################################# GEV quantile function
Q <- function(p, mu, sigma, xi) {
  if (xi == 0) {
    mu - sigma * log(-log(p))
  } else {
    mu + sigma / xi * (1/(-log(p))^xi - 1)
  }
}

##################################################
n0=function(y,varphis){
  nn=length(varphis)
  n0s=rep(NA,nn)
  for(i in 1:nn){
    n0s[i]=sum(abs(y)<=varphis[i])
  }
  return(n0s)
}


##################################################
groupL0Fdr<-function(trueindex,truenonindex,index.selected){
  return(length(intersect(truenonindex,index.selected))/max(1,length(index.selected)))
}

##################################################
groupL0Tpr<-function(trueindex,truenonindex,index.selected){
  if(length(trueindex)!=0){
    res=length(intersect(trueindex,index.selected))/length(trueindex)
  }else{
    res=1
  }
  return(res)
}

##################################################
groupL1Fdr<-function(trueindex,truenonindex,index.selected,alpha.true,alpha.hat){
  l1norm.alpha.hat=apply(abs(alpha.hat),2,sum)
  inter.set=intersect(truenonindex,index.selected)
  if(length(inter.set)==0){
    res=0
    }else{
      res=sum(l1norm.alpha.hat[inter.set])/sum(l1norm.alpha.hat)##note zero at den
    }
  
  return(res)
}

##################################################
groupL1Tpr<-function(trueindex,truenonindex,index.selected,alpha.true,alpha.hat){
  l1norm.alpha.true=apply(abs(alpha.true),2,sum)
  l1norm.alpha.hat=apply(abs(alpha.hat),2,sum)
  inter.set=intersect(trueindex,index.selected)
  if(length(inter.set)==0){
    res=0
    }else{
      res=sum(l1norm.alpha.hat[inter.set])/sum(l1norm.alpha.true)##note zero at den
    }
  return(res)
}

##################################################
L1Fdr<-function(alpha.true,alpha.hat){
  alpha.true.vec=as.vector(alpha.true)
  alpha.hat.vec=as.vector(alpha.hat)
  trueindex=which(alpha.true!=0)
  truenonindex=which(alpha.true==0)
  index.selected=which(alpha.hat.vec!=0)
  inter.set=intersect(truenonindex,index.selected)
  if(length(inter.set)==0){
    res=0
  }else{
    res=sum(abs(alpha.hat.vec)[inter.set])/sum(abs(alpha.hat.vec))######note zero at den
  }
  
  return(res)
}

##################################################
L1Tpr<-function(alpha.true,alpha.hat){
  alpha.true.vec=as.vector(alpha.true)
  alpha.hat.vec=as.vector(alpha.hat)
  trueindex=which(alpha.true!=0)
  truenonindex=which(alpha.true==0)
  index.selected=which(alpha.hat.vec!=0)
  inter.set=intersect(trueindex,index.selected)
  if(length(inter.set)==0){
    res=0
  }else{
    res=sum(abs(alpha.hat.vec)[inter.set])/sum(abs(alpha.true.vec))######note zero at den
  }
  
  return(res)
}

#################################################
softfunction<-function(x,lambda){
  if( x > lambda){
    out=x-lambda
  }else if(x<(-lambda)){
      out=x+lambda
  }else{
      out=0
  }
  return(out)
}






