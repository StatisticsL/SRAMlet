#############################################################################
######------------simulations for  p=1-------------------------
rm(list = ls())
setwd("  ")##add your directory
source('loadAll.R')

##########################################################################
####------------parameter settings----------------------------------------
Ntrain=2^10#####a power of 2
Ntest=Ntrain
P=1
p0=2^3###a power of 2, the number of father wavelets
J=log2(p0)
SD=1 ##standard deviation
SNR=3 ###the parameter in wave-functions
num_simulations=100
len.lambda.grids=100

family="DaubExPhase";
fun.type='waveblocks'##options:wavedoppler,waveblocks,wavebumps,waveheavisine
###waveblocks filter.number=1
###others filter.number=3
filter.number=1;

bc="periodic";
alpha=0.05;
M_MC=1000;
M_GEV=100
type='GEV'
pen.father=F

###########################################################################
####-------------results records--------------------------------------------
mse=matrix(NA,nrow = num_simulations,ncol = 6)
mse=as.data.frame(mse)
names(mse)=c('sqrtQUT','sqrtsure','sqrtmse','softQUT','softsure','softmse')
fdr=mse
tpr=mse

num_needles=rep(0,num_simulations)
trueindex_list=vector('list',length = num_simulations)

Lambdas.soft=vector('list',length = num_simulations)
Lambdas.sqrt=vector('list',length = num_simulations)

Mse.soft=vector('list',length = num_simulations)
Mse.sqrt=vector('list',length = num_simulations)

Sure.soft=vector('list',length = num_simulations)
Sure.sqrt=vector('list',length = num_simulations)

fitResults.sqrt=vector('list',length = num_simulations)
fitResults.soft=vector('list',length = num_simulations)

sigma.all=matrix(NA,nrow = num_simulations,ncol = 6)
sigma.all=as.data.frame(sigma.all)
names(sigma.all)=c('sqrtmse','sqrtQUT','sqrtsure','softmse','softQUT','softsure')
###################################################
for (simu in 1:num_simulations) {
  print(paste(simu,'begins',sep = ' '))
  ###########################################################################
  ####-------------generate data---------------------------------------------
  x.train=matrix(runif(Ntrain*P),Ntrain,P)
  x.test=matrix(runif(Ntest*P),Ntest,P)
  
  if(fun.type=='waveheavisine'){
    y.train=waveheavisine(x=x.train,snr = SNR)
    y.test=waveheavisine(x=x.test,snr = SNR)
  }else if(fun.type=='wavedoppler'){
    y.train=wavedoppler(x=x.train,snr = SNR)
    y.test=wavedoppler(x=x.test,snr = SNR)
  }else if(fun.type=='waveblocks'){
    y.train=waveblocks(x=x.train,snr = SNR)
    y.test=waveblocks(x=x.test,snr = SNR)
  }else if(fun.type=='wavebumps'){
    y.train=wavebumps(x=x.train,snr = SNR)
    y.test=wavebumps(x=x.test,snr = SNR)
  }
  mu.train=y.train
  y.train=y.train+rnorm(Ntrain,sd=SD)
  
  alpha.true=rep(0,Ntrain)
  mut=wd(mu.train[order(x.train)],filter.number=filter.number, family=family, bc=bc)
  alpha.true[1:p0]=accessC(mut,J)
  ii=p0
  for(level in J:(mut$nlevels-1)){
    lj=2^level
    alpha.true[(ii+1):(ii+lj)]=accessD(mut,level)
    ii=ii+lj
  }
  
  trueindex=setdiff(which(alpha.true !=0),seq(1:p0))
  truenonindex=setdiff(which(alpha.true ==0),seq(1:p0))
  num_needles[simu]=length(trueindex)
  trueindex_list[[simu]]= trueindex
  
  # trueindex=which(abs(alpha.true)/sum(abs(alpha.true)) > 0.001)
  # truenonindex=which(abs(alpha.true)/sum(abs(alpha.true)) <= 0.001)
  
  # filename.data=paste('./simulation1/MSE/raw data/data-',fun.type,'-N-',Ntrain,'-P0-',p0,'-SNR-',SNR,'-simu-',simu,'.RData',sep = '')
  # save(x.train,y.train,x.test,y.test,mu.train,file = filename.data)  
  ############################################################################
  ####-------------QUT----------------------
  fit_soft.QUT=softwavethresh(y=y.train,X=x.train,lambda=NA,sigma=NA,p0=p0,pen.father=pen.father,
                              filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                              type=type,com.sure=T,com.sqrtrss=F)
  fit_sqrt.QUT=sqrtwavethresh(y=y.train,X=x.train,lambda=NA,sigma=NA, p0=p0,pen.father=pen.father,
                              filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                              type=type,com.sure=T,com.sqrtrss=T)
  
  sigma.all[simu,'softQUT']=fit_soft.QUT$sigma
  sigma.all[simu,'sqrtQUT']=fit_sqrt.QUT$sqrtrss
  
  selectedindex.sqrt.QUT=which(abs(fit_sqrt.QUT$alpha.hat) !=0)
  selectedindex.soft.QUT=which(abs(fit_soft.QUT$alpha.hat) !=0)
  
  fdr[simu,'sqrtQUT']=length(intersect(truenonindex,selectedindex.sqrt.QUT))/max(length(selectedindex.sqrt.QUT),1)
  fdr[simu,'softQUT']=length(intersect(truenonindex,selectedindex.soft.QUT))/max(length(selectedindex.soft.QUT),1)
  tpr[simu,'sqrtQUT']=length(intersect(trueindex,selectedindex.sqrt.QUT))/length(trueindex)
  tpr[simu,'softQUT']=length(intersect(trueindex,selectedindex.soft.QUT))/length(trueindex)
  
  
  lambdaQUT.soft=fit_soft.QUT$lambdaQUT
  lambdaQUT.sqrt=fit_sqrt.QUT$lambdaQUT
  
  y.test.hat.QUT.soft=predict.softwavethresh(x.train,fit_soft.QUT,x.test)
  y.test.hat.QUT.sqrt=predict.sqrtwavethresh(x.train,fit_sqrt.QUT,x.test)
  
  mse[simu,'sqrtQUT']=mean((y.test-y.test.hat.QUT.sqrt)^2)
  mse[simu,'softQUT']=mean((y.test-y.test.hat.QUT.soft)^2)
  ############################################################################
  ####-------------MSE and SURE----------------------
  Lambdas.soft[[simu]]=exp(seq(log(lambdaQUT.soft/10),log(lambdaQUT.soft*1.2),length=len.lambda.grids))
  Lambdas.sqrt[[simu]]=exp(seq(log(lambdaQUT.sqrt/10),log(lambdaQUT.sqrt*1.2),length=len.lambda.grids))
  
  Mse.soft[[simu]]=rep(NA,length(Lambdas.soft[[simu]]))
  Mse.sqrt[[simu]]=rep(NA,length(Lambdas.sqrt[[simu]]))
  
  Sure.soft[[simu]]=rep(NA,length(Lambdas.soft[[simu]]))
  Sure.sqrt[[simu]]=rep(NA,length(Lambdas.sqrt[[simu]]))
  
  for (i in 1:len.lambda.grids) {
    fit_soft.temp=softwavethresh(y=y.train,X=x.train,lambda=Lambdas.soft[[simu]][i],sigma=NA,p0=p0,pen.father=pen.father,
                                 filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                                 type=type,com.sure=T,com.sqrtrss=F)
    fit_sqrt.temp=sqrtwavethresh(y=y.train,X=x.train,lambda=Lambdas.sqrt[[simu]][i],sigma=NA, p0=p0,pen.father=pen.father,
                                 filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                                 type=type,com.sure=T,com.sqrtrss=F) 
    y.train.hat.temp.soft=predict.softwavethresh(x.train,fit_soft.temp,x.train)
    Mse.soft[[simu]][i]=sum((mu.train-y.train.hat.temp.soft)^2)
    y.train.hat.temp.sqrt=predict.sqrtwavethresh(x.train,fit_sqrt.temp,x.train)
    Mse.sqrt[[simu]][i]=sum((mu.train-y.train.hat.temp.sqrt)^2)
    
    Sure.soft[[simu]][i]=fit_soft.temp$SURE
    Sure.sqrt[[simu]][i]=fit_sqrt.temp$SURE
  }
  
  lambda.mse.soft=Lambdas.soft[[simu]][which.min(Mse.soft[[simu]])]
  lambda.mse.sqrt=Lambdas.sqrt[[simu]][which.min(Mse.sqrt[[simu]])]
  
  lambda.sure.soft=Lambdas.soft[[simu]][which.min(Sure.soft[[simu]])]
  lambda.sure.sqrt=Lambdas.sqrt[[simu]][which.min(Sure.sqrt[[simu]])]
  
  fit_soft.mse=softwavethresh(y=y.train,X=x.train,lambda=lambda.mse.soft,sigma=NA,p0=p0,pen.father=pen.father,
                              filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                              type=type,com.sure=T,com.sqrtrss=F)
  fit_sqrt.mse=sqrtwavethresh(y=y.train,X=x.train,lambda= lambda.mse.sqrt,sigma=NA, p0=p0,pen.father=pen.father,
                              filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                              type=type,com.sure=T,com.sqrtrss=F) 
  fit_soft.sure= softwavethresh(y=y.train,X=x.train,lambda=lambda.sure.soft,sigma=NA,p0=p0,pen.father=pen.father,
                                filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                                type=type,com.sure=T,com.sqrtrss=F)
  fit_sqrt.sure=sqrtwavethresh(y=y.train,X=x.train,lambda=lambda.sure.sqrt,sigma=NA, p0=p0,pen.father=pen.father,
                               filter.number=filter.number,family=family, bc=bc,alpha=alpha,M_MC=M_MC,M_GEV=M_GEV,
                               type=type,com.sure=T,com.sqrtrss=F)
  y.test.hat.mse.soft=predict.softwavethresh(x.train,fit_soft.mse,x.test)
  y.test.hat.mse.sqrt=predict.sqrtwavethresh(x.train,fit_sqrt.mse,x.test)
  y.test.hat.sure.soft=predict.softwavethresh(x.train,fit_soft.sure,x.test)
  y.test.hat.sure.sqrt=predict.sqrtwavethresh(x.train,fit_sqrt.sure,x.test)
  
  mse[simu,'sqrtmse']=mean((y.test-y.test.hat.mse.sqrt)^2)
  mse[simu,'sqrtsure']=mean((y.test-y.test.hat.sure.sqrt)^2)
  mse[simu,'softmse']=mean((y.test-y.test.hat.mse.soft)^2)
  mse[simu,'softsure']=mean((y.test-y.test.hat.sure.soft)^2)
  
  selectedindex.sqrt.mse=setdiff(which(abs(fit_sqrt.mse$alpha.hat) !=0),seq(1,p0))
  selectedindex.soft.mse=setdiff(which(abs(fit_soft.mse$alpha.hat) !=0),seq(1,p0))
  fdr[simu,'sqrtmse']=length(intersect(truenonindex,selectedindex.sqrt.mse))/max(length(selectedindex.sqrt.mse),1)
  fdr[simu,'softmse']=length(intersect(truenonindex,selectedindex.soft.mse))/max(length(selectedindex.soft.mse),1)
  tpr[simu,'sqrtmse']=length(intersect(trueindex,selectedindex.sqrt.mse))/length(trueindex)
  tpr[simu,'softmse']=length(intersect(trueindex,selectedindex.soft.mse))/length(trueindex)
  
  selectedindex.sqrt.sure=setdiff(which(abs(fit_sqrt.sure$alpha.hat) !=0),seq(1,p0))
  selectedindex.soft.sure=setdiff(which(abs(fit_soft.sure$alpha.hat) !=0),seq(1,p0))
  fdr[simu,'sqrtsure']=length(intersect(truenonindex,selectedindex.sqrt.sure))/max(length(selectedindex.sqrt.sure),1)
  fdr[simu,'softsure']=length(intersect(truenonindex,selectedindex.soft.sure))/max(length(selectedindex.soft.sure),1)
  tpr[simu,'sqrtsure']=length(intersect(trueindex,selectedindex.sqrt.sure))/length(trueindex)
  tpr[simu,'softsure']=length(intersect(trueindex,selectedindex.soft.sure))/length(trueindex)
  
  fitResults.sqrt[[simu]][1]= fit_sqrt.QUT
  fitResults.sqrt[[simu]][2]= fit_sqrt.mse
  fitResults.sqrt[[simu]][3]= fit_sqrt.sure
  
  fitResults.soft[[simu]][1]= fit_soft.QUT
  fitResults.soft[[simu]][2]= fit_soft.mse
  fitResults.soft[[simu]][3]= fit_soft.sure
}

# ###############################################################################
# ####-------------save results--------------------------------------------
# save(mse,Lambdas.soft,Lambdas.sqrt,Mse.soft,Mse.sqrt,Sure.soft,Sure.sqrt,
#      file = paste('./simulation1/MSE/results/compution-',fun.type,'-N-',Ntrain,'-P0-',p0,'-SNR-',SNR,'.RData',sep = ''))
# save(fitResults.sqrt,fitResults.soft,
#      file =paste('./simulation1/MSE/results/fitting-',fun.type,'-N-',Ntrain,'-P0-',p0,'-SNR-',SNR,'.RData',sep = ''))
# 

##########################################################################
####------------------first plot----------------
##############################################################
par(mfrow=c(2,2))

###here is the oracle mse
y.train.hat.QUT.soft=predict.softwavethresh(x.train,fit_soft.QUT,x.train)
y.train.hat.QUT.sqrt=predict.sqrtwavethresh(x.train,fit_sqrt.QUT,x.train)
Mse.QUT.soft=sum((mu.train-y.train.hat.QUT.soft)^2)
Mse.QUT.sqrt=sum((mu.train-y.train.hat.QUT.sqrt)^2)
Sure.QUT.soft=fit_soft.QUT$SURE
Sure.QUT.sqrt=fit_sqrt.QUT$SURE

lambdas1=Lambdas.soft[[num_simulations]]/lambdaQUT.soft
lambdas2=Lambdas.sqrt[[num_simulations]]/lambdaQUT.sqrt


#######################################################################
mainname=paste('Square-root','       ','Original',sep = '')

boxplot(fdr,main=mainname,par(las="2"),xaxt='n',ylab='FDR',ylim=c(0,1))##par(las="2")坐标竖着显示
axis(1,at=c(1,2,3,4,5,6),labels = c('QUT','SURE','oracle','QUT','SURE','oracle'))
abline(v=3.5,lwd=1)

boxplot(tpr,main=mainname,par(las="2"),xaxt='n',ylab='TPR',ylim=c(0,1))##par(las="2")坐标竖着显示
axis(1,at=c(1,2,3,4,5,6),labels = c('QUT','SURE','oracle','QUT','SURE','oracle'))
abline(v=3.5,lwd=1)

boxplot(mse,main=mainname,par(las="2"),xaxt='n',ylab=expression(l[2]-loss))##par(las="2")坐标竖着显示
axis(1,at=c(1,2,3,4,5,6),labels =c('QUT','SURE','oracle','QUT','SURE','oracle'))
abline(v=3.5,lwd=1)

boxplot(sigma.all,par(las="1"),xaxt='n',main=mainname,ylab=expression(hat(bold(sigma))))#
#axis(1,at=c(2,5),labels = c('QUT','QUT'))# c(expression(Square-root[QUT]),expression(original[QUT]))
abline(v=3.5,lwd=1)
abline(h=1,lty=4,lwd=1)



####----------------second plot-------------------------------------
par(mfrow=c(2,2),mai=c(0.8,0.8,0.8,0.2))
y.max=max(Mse.sqrt[[num_simulations]]/Ntrain,Sure.sqrt[[num_simulations]]/Ntrain,Mse.soft[[num_simulations]]/Ntrain,Sure.soft[[num_simulations]]/Ntrain)
plot(lambdas2,Mse.sqrt[[num_simulations]]/Ntrain,type='l',xlim=c(min(lambdas1,lambdas2),max(lambdas1,lambdas2)),
     ylim = c(0,y.max),xlab = expression(lambda/lambda[QUT]),
     ylab = expression(paste(l[2]-loss,' and SURE',sep = ' ')),main = expression(Blocks))
points(1,Mse.QUT.sqrt/Ntrain,pch=16)
lines(lambdas2,Sure.sqrt[[num_simulations]]/Ntrain,lty=2)
# points(1,Sure.QUT.sqrt/Ntrain,pch=16)
lines(lambdas1,Mse.soft[[num_simulations]]/Ntrain,col='2')
points(1,Mse.QUT.soft/Ntrain,pch=16,col='2')
lines(lambdas1,Sure.soft[[num_simulations]]/Ntrain,lty=2,col='2')
# points(1,Sure.QUT.soft/Ntrain,pch=16,col='2')


plot(x.test[,1][order(x.test[,1])],y.test[,1][order(x.test[,1])],type = 'l',col='blue',main = expression(l[2]-loss),xlab = 'x',ylab = expression(hat(mu)))
lines(x.test[,1][order(x.test[,1])],y.test.hat.mse.soft[order(x.test[,1])],col='red')
lines(x.test[,1][order(x.test[,1])],y.test.hat.mse.sqrt[order(x.test[,1])])


plot(x.test[,1][order(x.test[,1])],y.test[,1][order(x.test[,1])],type = 'l',col='blue',main = expression(SURE),xlab = 'x',ylab = expression(hat(mu)))
lines(x.test[,1][order(x.test[,1])],y.test.hat.sure.soft[order(x.test[,1])],col='red')
lines(x.test[,1][order(x.test[,1])],y.test.hat.sure.sqrt[order(x.test[,1])])

plot(x.test[,1][order(x.test[,1])],y.test[,1][order(x.test[,1])],type = 'l',col='blue',main = expression(QUT),xlab = 'x',ylab = expression(hat(mu)))
lines(x.test[,1][order(x.test[,1])],y.test.hat.QUT.soft[order(x.test[,1])],col='red')
lines(x.test[,1][order(x.test[,1])],y.test.hat.QUT.sqrt[order(x.test[,1])])




