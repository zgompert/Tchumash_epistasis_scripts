library(BMS)
load("models.rdat")

## leave-one-out cross-validation
pc<-prcomp(c_AC_pip[,-6],center=TRUE,scale=FALSE)

pc_AC<-pc$x[,1]

V<-15*2 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
        pcmod[,i]<-pcmod[,i] * pc_AC
}
dfAC<-data.frame(fit_AC$exfit,blk2=as.numeric(fit_AC$blk==2),blk3=as.numeric(fit_AC$blk==3),mod,pcmod)
Ni<-dim(dfAC)[1]
kgrp<-sample(1:Ni,Ni,replace=FALSE)
predAC<-rep(NA,Ni)
for(i in 1:200){
        dfACx<-dfAC[-c(kgrp==i),]
        o<-bms(dfACx,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
        j<-which(kgrp==i)
        predAC[j]<-sum(coef(o,order.by.pip=FALSE)[,2]*dfAC[j,-1])
}
cor.test(dfAC[,1],predAC)

pc<-prcomp(c_MM_pip[,-6],center=TRUE,scale=FALSE)
pc_MM<-pc$x[,1]

form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
        pcmod[,i]<-pcmod[,i] * pc_MM
}

dfMM<-data.frame(fit_MM$exfit,blk2=as.numeric(fit_MM$blk==2),blk3=as.numeric(fit_MM$blk==3),mod,pcmod)
Ni<-dim(dfMM)[1]
kgrp<-sample(1:Ni,Ni,replace=FALSE)
predMM<-rep(NA,Ni)
for(i in 1:200){
        dfMMx<-dfMM[-c(kgrp==i),]
        o<-bms(dfMMx,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
        j<-which(kgrp==i)
        predMM[j]<-sum(coef(o,order.by.pip=FALSE)[,2]*dfMM[j,-1])
}
cor.test(dfMM[,1],predMM)

save(list=ls(),file="models_cv.rdat")

## remove block effect
mnObsAC<-tapply(INDEX=fit_AC$blk,X=dfAC[,1],mean,na.rm=TRUE)
mnPredAC<-tapply(INDEX=fit_AC$blk,X=predAC,mean,na.rm=TRUE)
rObsAC<-dfAC[,1]
rPredAC<-predAC
for(i in 1:3){
	a<-which(fit_AC$blk==i)
	rObsAC[a]<-rObsAC[a]-mnObsAC[i]
	rPredAC[a]<-rPredAC[a]-mnPredAC[i]
}

cor.test(rObsAC,rPredAC)
#t = 9.567, df = 195, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.4618854 0.6533604
#sample estimates:
#      cor 
#0.5651872 


mnObsMM<-tapply(INDEX=fit_MM$blk,X=dfMM[,1],mean,na.rm=TRUE)
mnPredMM<-tapply(INDEX=fit_MM$blk,X=predMM,mean,na.rm=TRUE)
rObsMM<-dfMM[,1]
rPredMM<-predMM
for(i in 1:3){
	a<-which(fit_MM$blk==i)
	rObsMM[a]<-rObsMM[a]-mnObsMM[i]
	rPredMM[a]<-rPredMM[a]-mnPredMM[i]
}

cor.test(rObsMM,rPredMM)

#t = 7.9074, df = 195, p-value = 1.907e-13
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.3790579 0.5917751
#sample estimates:
#      cor 
#0.4927429 
