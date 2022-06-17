rm(list = ls())
setwd("  ")##add your directory

MSE=matrix(NA,nrow = 100*4*2,ncol=2)
MSE=as.data.frame(MSE)
names(MSE)=c('Values','Methods','N')
FDR=MSE
TPR=MSE

J.vec=c(8,9,10,11)
sigma.vec.SR=rep(0,4)
sigma.vec.MAD=rep(0,4)

SNR=3
S=4

for (jj in 1:4) {
  
  j=jj+7
  N=2^(j)
  P=2^(j+1)
  load(paste('results/compution-N-',N,'-P-',P,'-S-',S,'-SNR-',SNR,'.RData',sep = ''))
  sigma.vec.SR[jj]=mean(sigma.all[,'sqrtAMlet'])
  sigma.vec.MAD[jj]=mean(sigma.all[,'madAMlet'])
  
  namemse.vec=colnames(mse)
  namefdr.vec=colnames(fdr.group.l0)
  nametpr.vec=colnames(tpr.group.l0)
  
  jn=(jj-1)*200
  MSE[(jn+1):(jn+200),'N']=N
  FDR[(jn+1):(jn+200),'N']=N
  TPR[(jn+1):(jn+200),'N']=N
  for (i in 1:2) {
    MSE[(jn+(i-1)*100+1):(jn+(i-1)*100+100),'Values']=mse[,namemse.vec[i]]
    MSE[(jn+(i-1)*100+1):(jn+(i-1)*100+100),'Methods']= namemse.vec[i]
    
    FDR[(jn+(i-1)*100+1):(jn+(i-1)*100+100),'Values']=fdr.group.l0[,namefdr.vec[i]]
    FDR[(jn+(i-1)*100+1):(jn+(i-1)*100+100),'Methods']= namefdr.vec[i]
    
    TPR[(jn+(i-1)*100+1):(jn+(i-1)*100+100),'Values']=tpr.group.l0[,nametpr.vec[i]]
    TPR[(jn+(i-1)*100+1):(jn+(i-1)*100+100),'Methods']= nametpr.vec[i]
  }
}

######----------plot MSE--------------------------------------------
par(mfrow=c(2,2))

######----------plot FDR and TPR--------------------------------------------
sqrtFDR=rep(0,4)
madFDR=rep(0,4)
sqrtFDR[1]=mean(FDR[1:100,'Values'])
madFDR[1]=mean(FDR[101:200,'Values'])
sqrtFDR[2]=mean(FDR[201:300,'Values'])
madFDR[2]=mean(FDR[301:400,'Values'])
sqrtFDR[3]=mean(FDR[401:500,'Values'])
madFDR[3]=mean(FDR[501:600,'Values'])
sqrtFDR[4]=mean(FDR[601:700,'Values'])
madFDR[4]=mean(FDR[701:800,'Values'])
plot(J.vec,sqrtFDR,type='o',xaxt='n',ylim = c(0,1),ylab = 'FDR',xlab = expression(paste('n=',2^j,', ','p=',2^(j+1))))
axis(1,at=c(8,9,10,11),labels = c('j=8','j=9','j=10','j=11'))
#abline(h=1,lty=4,lwd=1)
lines(J.vec,madFDR,col='red',type='o')


# 
# 
# boxplot(FDR$Values~FDR$Methods+FDR$N,col=rep(col.vec,4),xaxt='n',ylab = 'FDR',xlab = expression(paste('n=',2^j,', ','p=',2^(j+1))))#,par(las="1"),
# axis(1,at=c(1.5,3.5,5.5,7.5),labels = c('j=8','j=9','j=10','j=11'))
# # letters=c('sqrtAMlet','oracleAMlet','madAMlet')
# # legend(x=c(9,12),y=c(0.8,1),legend=letters, cex=0.8,col=c('red','green','blue'),pch=16)  

sqrtTPR=rep(0,4)
madTPR=rep(0,4)
sqrtTPR[1]=mean(TPR[1:100,'Values'])
madTPR[1]=mean(TPR[101:200,'Values'])
sqrtTPR[2]=mean(TPR[201:300,'Values'])
madTPR[2]=mean(TPR[301:400,'Values'])
sqrtTPR[3]=mean(TPR[401:500,'Values'])
madTPR[3]=mean(TPR[501:600,'Values'])
sqrtTPR[4]=mean(TPR[601:700,'Values'])
madTPR[4]=mean(TPR[701:800,'Values'])
plot(J.vec,sqrtTPR,type='o',xaxt='n',ylim = c(0,1),ylab = 'TPR',xlab = expression(paste('n=',2^j,', ','p=',2^(j+1))))
axis(1,at=c(8,9,10,11),labels = c('j=8','j=9','j=10','j=11'))
#abline(h=1,lty=4,lwd=1)
lines(J.vec,madTPR,col='red',type = 'o')


# boxplot(TPR$Values~TPR$Methods+TPR$N,col=rep(col.vec,4),xaxt='n',ylab = 'TPR',xlab = expression(paste('n=',2^j,', ','p=',2^(j+1))))#,par(las="1"),
# axis(1,at=c(1.5,3.5,5.5,7.5),labels = c('j=8','j=9','j=10','j=11'))
# # letters=c('sqrtAMlet','oracleAMlet','madAMlet')
# # legend(x=c(8,12.5),y=c(0.05,0.2),legend=letters, cex=0.8,col=c('red','green','blue'),pch=16)
# # 

# J=c(7,8,9,10,11,12,13,14,15)
# SR=c(5.953389,5.883509,5.480874,4.831618,3.927436,2.801428,1.925522,1.481096, 1.248860)
# MA=c(5.351196 ,5.169868,4.738971,4.344797,3.894915,3.480085,3.104179, 2.782671,2.530548)
# plot(J,SR,type='l',ylim = c(0,7),main='P=10',xaxt='n',xlab = 'n',ylab = expression(hat(sigma)))
# axis(1,at=c(8,10,12,14),labels = c(expression(2^8),expression(2^10),expression(2^12),expression(2^14)))
# axis(2,at=1,labels = 1)
# abline(h=1,lty=4,lwd=1)
# lines(J,MA,col='red')


col.vec=c('red','white')
boxplot(MSE$Values~MSE$Methods+MSE$N,xaxt='n',col=rep(col.vec,4),ylab = expression(l[2]-loss),xlab = expression(paste('n=',2^j,', ','p=',2^(j+1))))#,par(las="1"),
axis(1,at=c(1.5,3.5,5.5,7.5),labels = c('j=8','j=9','j=10','j=11'))
letters=c('sqrtAMlet','madAMlet')
legend(x=c(9,12),y=c(38,45),legend=letters, col=col.vec,cex=1,pch=16)


plot(J.vec,sigma.vec.SR,type='o',xaxt='n',ylim = c(0,7),ylab = expression(hat(sigma)),xlab = expression(paste('n=',2^j,', ','p=',2^(j+1))))
axis(1,at=c(8,9,10,11),labels = c('j=8','j=9','j=10','j=11'))
abline(h=1,lty=4,lwd=1)
lines(J.vec,sigma.vec.MAD,col='red',type='o')



