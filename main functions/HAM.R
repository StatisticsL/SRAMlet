HAM<-function(x.train,y.train,x.val,y.val,lambda.pen,lambda.curv,knots, model = LinReg(),control = grpl.control(trace = 0)){
  fit_HAM=penGAM(x=x.train,y=y.train,lambda.pen =lambda.pen,
                 lambda.curv = lambda.curv ,knots = knots,
                 model = LinReg(),control = grpl.control(trace = 0))
  pred_HAM=predict(fit_HAM,newdata = x.val,type = 'response')
  
  pen.length=length(lambda.pen)
  curv.length=length(lambda.curv)
  
  # mse.vec=rep(0,pen.length*curv.length)
  # mse.fun<-function(x,y){
  #   return(mean((x-y)^2))
  # }
  # for (i in 1:pen.length) {
  #   mse.vec[((i-1)*curv.length+1):((i-1)*curv.length+curv.length)]=apply(pred_HAM[i,,],1,mse.fun,y=y.val)
  # }
  # mse.Matrix=t(matrix(mse.vec,nrow=curv.length,ncol=pen.length))
  
  mse.Matrix=matrix(0,pen.length,curv.length)
  for(i in 1:pen.length){
    for (j in 1:curv.length) {
      mse.Matrix[i,j]=mean((pred_HAM[i,j,]-y.val)^2)
    }
  }
  index=which(mse.Matrix==min(mse.Matrix),arr.ind = TRUE)
  index.pen=index[1]
  index.curv=index[2]
  
  coeff=fit_HAM$coefficients[index.pen,index.curv,]
  coeff.Matrix=matrix(coeff[-1],nrow=(knots+2))
  index.selected=which(colSums(abs(coeff.Matrix))!=0)
  
  out=NULL
  out$coeff.Matrix=coeff.Matrix
  out$mse.Matrix=mse.Matrix
  out$index.pen=index.pen
  out$index.curv=index.curv
  out$index.selected=index.selected
  out$pred.val=pred_HAM
  return(out)
}