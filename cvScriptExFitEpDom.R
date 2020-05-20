library(BMS)
load("models.rdat")

## leave-one-out cross-validation

V<-33 ## variables
form<-p_AC[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
d3<-1-abs(g_AC_pip[,3]-1)## het = 1, both homos = 0
d4<-1-abs(g_AC_pip[,4]-1)

dfAC<-data.frame(fit_AC$exfit,blk2=as.numeric(fit_AC$blk==2),blk3=as.numeric(fit_AC$blk==3),mod,d3,d4)
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

form<-p_MM[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
d3<-1-abs(g_MM_pip[,3]-1)## het = 1, both homos = 0
d4<-1-abs(g_MM_pip[,4]-1)

dfMM<-data.frame(fit_MM$exfit,blk2=as.numeric(fit_MM$blk==2),blk3=as.numeric(fit_MM$blk==3),mod,d3,d4)
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

save(list=ls(),file="models_cvEpiDom.rdat")

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
