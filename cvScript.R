library(BMS)
load("models.rdat")

## leave-one-out cross-validation
V<-15*3 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/3)){
        pcmod[,i]<-pcmod[,i] * pc_AC
}
msmod<-mod
for(i in 1:(V/3)){
        msmod[,i]<-msmod[,i] * ms_AC
}
dfAC<-data.frame(p_AC[,2],mod,pcmod,msmod)

Ni<-dim(dfAC)[1]
kgrp<-sample(1:Ni,Ni,replace=FALSE)
predAC<-rep(NA,Ni)
for(i in 1:200){
        dfACx<-dfAC[-c(kgrp==i),]
        o<-bms(dfACx,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
        j<-which(kgrp==i)
        predAC[j]<-sum(coef(o,order.by.pip=FALSE)[,2]*dfAC[j,-1])
}
cor.test(dfAC[,1],predAC)

form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/3)){
        pcmod[,i]<-pcmod[,i] * pc_MM
}
msmod<-mod
for(i in 1:(V/3)){
        msmod[,i]<-msmod[,i] * ms_MM
}
dfMM<-data.frame(p_MM[,2],mod,pcmod,msmod)

Ni<-dim(dfMM)[1]
kgrp<-sample(1:Ni,Ni,replace=FALSE)
predMM<-rep(NA,Ni)
for(i in 1:200){
        dfMMx<-dfMM[-c(kgrp==i),]
        o<-bms(dfMMx,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
        j<-which(kgrp==i)
        predMM[j]<-sum(coef(o,order.by.pip=FALSE)[,2]*dfMM[j,-1])
}
cor.test(dfMM[,1],predMM)
cov(dfMM[,1],predMM)/var(dfMM[,1])

save(list=ls(),file="models.rdat")
