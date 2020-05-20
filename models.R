## read in genetic and survival data
g_AC<-as.matrix(read.table("g_tchum_AC.txt",header=FALSE))
g_MM<-as.matrix(read.table("g_tchum_MM.txt",header=FALSE))

snps<-as.matrix(read.table("tchumSnpTable.txt",header=FALSE))
del<-which(snps[,2]==128 & snps[,3] >= 5000000 & snps[,3] <= 6000000)

g_AC_del<-t(g_AC[del,])
g_MM_del<-t(g_MM[del,])


c_AC_del<-g_AC_del
c_MM_del<-g_MM_del
for(i in 1:length(del)){
	c_AC_del[,i]<-c_AC_del[,i]-mean(c_AC_del[,i])
	c_MM_del[,i]<-c_MM_del[,i]-mean(c_MM_del[,i])
}

c_AC_del2<-c_AC_del^2
c_MM_del2<-c_MM_del^2

ms<-c(which(snps[,2]==128 & snps[,3] <5000000),which(snps[,2]==702.1 & snps[,3] < 4139489))
g_AC_ms<-t(g_AC[ms,])
g_MM_ms<-t(g_MM[ms,])

c_AC_ms<-g_AC_ms
c_MM_ms<-g_MM_ms
for(i in 1:length(ms)){
	c_AC_ms[,i]<-c_AC_ms[,i]-mean(c_AC_ms[,i])
	c_MM_ms[,i]<-c_MM_ms[,i]-mean(c_MM_ms[,i])
}

p_AC<-read.table("pheno_AC.txt",header=FALSE)
p_MM<-read.table("pheno_MM.txt",header=FALSE)

## top pip SNPs

f<-list.files("output/",pattern="pip_")
pips<-matrix(NA,nrow=11855,ncol=6)
for(i in 1:6){
        pips[,i]<-scan(paste("output/",pf[i],sep=""))
}

pipSnpTab<-read.table("output/pipSnpTable.txt",header=FALSE)
pipDel<-which(pipSnpTab[,2]==128 & pipSnpTab[,3] >= 5000000 & pipSnpTab[,3] <= 6000000)

pipDf<-data.frame(pipSnpTab[pipDel,],pips[pipDel,5:6])

sum(pipDf[,5] > 0.5 | pipDf[,4] > 0.5)
#[1] 6
cor.test(pipDf[,4],pipDf[,5])

#	Pearson's product-moment correlation

#Pearson's product-moment correlation

#data:  pipDf[, 4] and pipDf[, 5]
#t = 2.0706, df = 16, p-value = 0.05493
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.009121002  0.762851117
#sample estimates:
#      cor 
#0.4597067 
#data:  pipDf[, 4] and pipDf[, 5]

## get relevant SNP ids
topIds<-pipDf[which(pipDf[,5] > 0.5 | pipDf[,4] > 0.5),1:3]
topSnps<-rep(NA,6)
#
for(i in 1:6){
	topSnps[i]<-which(snps[,1]==topIds[i,1] & snps[,2]==topIds[i,2] & snps[,3]==topIds[i,3])
}

## RG = 2, 3, 5, 6
## GB = 1, 3, 4, 6

g_AC_pip<-t(g_AC[topSnps,])
g_MM_pip<-t(g_MM[topSnps,])

## get effect estimates
ef<-list.files("output/",pattern="mav_")
mav<-matrix(NA,nrow=11855,ncol=6)
for(i in 1:6){
        mav[,i]<-scan(paste("output/",ef[i],sep=""))
}


pips[as.numeric(row.names(topIds)),c(5:6)]
#         [,1]     [,2]
#[1,] 0.000237 0.712051
#[2,] 0.944754 0.013838
#[3,] 1.000000 0.877689
#[4,] 0.153342 0.729562
#[5,] 0.994310 0.000577
#[6,] 1.000000 1.000000
 mav[as.numeric(row.names(topIds)),c(5:6)]
#              [,1]          [,2]
#[1,]  1.023523e-05  2.136363e-02
#[2,] -2.978231e-02  4.242169e-04
#[3,] -6.836257e-02  2.591891e-02
#[4,] -3.198406e-03 -1.893820e-02
#[5,]  4.027592e-02 -7.155026e-05
#[6,] -6.979418e-02  1.055617e-01

# 1 = GB+, 2= RG-, 3=GB-RG+, 4=GB-, 5=RG+, 6=GB-RG+
# given param. estimates, assume GB- and RG+ are 'similar'
# thus flip genotypes 1 and 2
# then counting alleles that increase RG and decrease GB
g_AC_pip[,1]<-2-g_AC_pip[,1]
g_AC_pip[,2]<-2-g_AC_pip[,2]
g_MM_pip[,1]<-2-g_MM_pip[,1]
g_MM_pip[,2]<-2-g_MM_pip[,2]

c_AC_pip<-g_AC_pip
c_MM_pip<-g_MM_pip
c_AC_pip<-c_AC_pip-1
c_MM_pip<-c_MM_pip-1
c_AC_mean<-apply(c_AC_pip,2,mean)
c_MM_mean<-apply(c_MM_pip,2,mean)
# uncomment for actual centering
for(i in 1:length(topSnps)){
        c_AC_pip[,i]<-c_AC_pip[,i]-mean(c_AC_pip[,i])
        c_MM_pip[,i]<-c_MM_pip[,i]-mean(c_MM_pip[,i])
}

library(BMS)
V<-21 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip))
mod<-cbind(mod[,-1])
dfAC<-data.frame(p_AC[,2],mod)
bmsAC<-vector("list",3)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
	bpipsAC[,i]<-oo[,1]
	bmavAC[,i]<-oo[,2]
}

paIds<-row.names(oo)

form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip))
mod<-cbind(mod[,-1])
dfMM<-data.frame(p_MM[,2],mod)
bmsMM<-vector("list",3)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
	bpipsMM[,i]<-oo[,1]
	bmavMM[,i]<-oo[,2]
}

cs<-rep(c("cornsilk3","coral3"),c(6,15))

pdf("ColorSurvPips.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

pdf("ColorSurvModAvEffs.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.8,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.8,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

## example interaction SNP 1 x 3
est<-apply(bmavAC,1,mean)
b<-est[c(1,3,8)]
b# 8 is interaction
#-0.0035672881  0.0002590297  0.1774696682
cg1<-c(-1,0,1)-c_AC_mean[1]
cg2<-c(-1,0,1)-c_AC_mean[3]
gfit<-matrix(NA,nrow=3,ncol=3)
for(i in 1:3){for(j in 1:3){
	gfit[i,j]<-b[1]*cg1[i] + b[2]*cg2[j] + b[3]*cg1[i]*cg2[j]
}}

gfit
#            [,1]         [,2]        [,3]
#[1,]  0.37557553  0.066168499 -0.24323853
#[2,]  0.16015593  0.028218562 -0.10371880
#[3,] -0.05526368 -0.009731376  0.03580093



library(RColorBrewer)
cl<-1.3
cm<-1.3
cn<-1.3
cs<-brewer.pal(n=3,"BrBG")
pdf("epiExample.pdf",width=6,height=9)
par(mfrow=c(2,1))
barplot(gfit,beside=TRUE,xlab="Genotype",col=cs,names.arg=c("aa","Aa","AA"),ylab="Effect",cex.lab=cl,cex.names=cn)
title(main="Focal SNP 1",cex.main=cm)
legend(10.5,0.4,c("aa","Aa","AA"),fill=cs,bty='n')
#box()
barplot(t(gfit),beside=TRUE,xlab="Genotype",col=cs,names.arg=c("aa","Aa","AA"),ylab="Effect",cex.lab=cl,cex.names=cn)
title(main="Focal SNP 3",cex.main=cm)
#box()
dev.off()


##### version without snp 6


library(BMS)
V<-15 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
dfAC<-data.frame(p_AC[,2],mod)
bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
	bpipsAC[,i]<-oo[,1]
	bmavAC[,i]<-oo[,2]
}

paIds<-row.names(oo)

form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
dfMM<-data.frame(p_MM[,2],mod)
bmsMM<-vector("list",3)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
	bpipsMM[,i]<-oo[,1]
	bmavMM[,i]<-oo[,2]
}

cs<-rep(c("cornsilk3","coral3"),c(5,10))

#pdf("ColorSurvPips5snp.pdf",width=7,height=9)
pdf("ColorSurvPips5snpHet.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

pdf("ColorSurvModAvEffs5snpHet.pdf",width=7,height=9)
#pdf("ColorSurvModAvEffs5snp.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.2,.2),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.2,.2),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

## example interaction SNP 1 x 3
est<-apply(bmavAC,1,mean)
b<-est[c(1,3,8)]
b# 8 is interaction
#-0.0035672881  0.0002590297  0.1774696682
cg1<-c(-1,0,1)-c_AC_mean[1]
cg2<-c(-1,0,1)-c_AC_mean[3]
gfit<-matrix(NA,nrow=3,ncol=3)
for(i in 1:3){for(j in 1:3){
	gfit[i,j]<-b[1]*cg1[i] + b[2]*cg2[j] + b[3]*cg1[i]*cg2[j]
}}

gfit
#            [,1]         [,2]        [,3]
#[1,]  0.37557553  0.066168499 -0.24323853
#[2,]  0.16015593  0.028218562 -0.10371880
#[3,] -0.05526368 -0.009731376  0.03580093

##### version without snp 6 but with PC1
pc<-prcomp(c_AC_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.9763 0.4129 0.37197 0.30890 0.26942
#Proportion of Variance 0.6665 0.1192 0.09676 0.06673 0.05076
#Cumulative Proportion  0.6665 0.7857 0.88251 0.94924 1.00000

pc_AC<-pc$x[,1]

library(BMS)
V<-15*2 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
	pcmod[,i]<-pcmod[,i] * pc_AC
}
dfAC<-data.frame(p_AC[,2],mod,pcmod)
bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
	bpipsAC[,i]<-oo[,1]
	bmavAC[,i]<-oo[,2]
}

paIds<-row.names(oo)

pc<-prcomp(c_MM_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.911 0.4188 0.36091 0.34237 0.24600
#Proportion of Variance 0.632 0.1335 0.09918 0.08925 0.04608
#Cumulative Proportion  0.632 0.7655 0.86467 0.95392 1.00000

pc_MM<-pc$x[,1]

form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
	pcmod[,i]<-pcmod[,i] * pc_MM
}
dfMM<-data.frame(p_MM[,2],mod,pcmod)
bmsMM<-vector("list",5)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
	bpipsMM[,i]<-oo[,1]
	bmavMM[,i]<-oo[,2]
}

cs<-rep(c("cornsilk3","coral3","gold1","dodgerblue4"),c(5,10,5,10))

pdf("ColorSurvPips5snpPC.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

pdf("ColorSurvModAvEffs5snpPC.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.2,.2),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.2,.2),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

## number of additive effects, number of interactions in the model
et<-rep(c(1:4),c(5,10,5,10))
Nx<-10000
noAC<-matrix(NA,nrow=Nx,ncol=4)
noMM<-matrix(NA,nrow=Nx,ncol=4)
for(i in 1:Nx){
	sAC<-rbinom(n=30,size=1,prob=apply(bpipsAC,1,mean))
	sMM<-rbinom(n=30,size=1,prob=apply(bpipsMM,1,mean))
	noAC[i,]<-tapply(X=sAC,INDEX=et,sum)
	noMM[i,]<-tapply(X=sMM,INDEX=et,sum)
}

mnsd<-function(x=NA){
	c(mean(x),sd(x))
}

noQTL_AC<-matrix(NA,nrow=6,ncol=2)
noQTL_MM<-matrix(NA,nrow=6,ncol=2)
for(i in 1:4){
	noQTL_AC[i,]<-mnsd(noAC[,i])
	noQTL_MM[i,]<-mnsd(noMM[,i])
}
noQTL_AC[5,]<-mnsd(apply(noAC[,2:4],1,sum))
noQTL_MM[5,]<-mnsd(apply(noMM[,2:4],1,sum))
noQTL_AC[6,]<-mnsd(apply(noAC,1,sum))
noQTL_MM[6,]<-mnsd(apply(noMM,1,sum))

nms<-c("G","GxG","GxPC1","GxGxPC1","epis.","all")
pdf("ColorNoEffs.pdf",width=6,height=5)
par(mar=c(6,4.5,0.5,1))
cs<-c("orange","cadetblue")
ox<-barplot(t(as.matrix(cbind(noQTL_AC[,1],noQTL_MM[,1]))),ylim=c(0,6),col=cs,beside=TRUE,xlab="",ylab="No. of effects",cex.lab=1.5,las=2,cex.names=1.2)
noQTL<-cbind(noQTL_AC,noQTL_MM)
oxm<-apply(ox,2,mean)
axis(1,at=oxm,nms,las=2,cex.axis=1.2)
ox<-as.vector(ox)
ox<-ox[c(1,3,5,7,9,11,2,4,6,8,10,12)]
segments(ox,noQTL[,1]+noQTL[,2],ox,noQTL[,1]-noQTL[,2])
legend(1,5.8,c("AC","MM"),fill=cs,bty='n',cex=1.2)
box()
dev.off()

## example interaction SNP 1 x 2
est<-apply(bmavAC,1,mean)
b<-est[c(1,2,6,16,17,21)]
#-9.927926e-05 -1.577068e-02  2.697549e-02 -3.723190e-02 -3.850671e-03 6.150750e-02

cg1<-c(-1,0,1)-c_AC_mean[1]
cg2<-c(-1,0,1)-c_AC_mean[3]
pcs<-c(-1,0,1)

gfit<-vector("list",3)
for(k in 1:3){
	gfit[[k]]<-matrix(NA,nrow=3,ncol=3)

	for(i in 1:3){for(j in 1:3){
		gfit[[k]][i,j]<-
		b[1]*cg1[i] + b[2]*cg2[j] + b[3]*cg1[i]*cg2[j]+
		b[4]*cg1[i]*pcs[k] + b[5]*cg2[j]*pcs[k]+
		b[6]*cg1[i]*cg2[j]*pcs[k]
	}}
}
gfit

#[[1]]
#            [,1]        [,2]         [,3]
#[1,] -0.12249157 -0.07415683 -0.025822082
#[2,] -0.04413678 -0.03033404 -0.016531308
#[3,]  0.03421802  0.01348874 -0.007240533

#[[2]]
#           [,1]        [,2]         [,3]
#[1,] 0.07518786 0.012347758 -0.050492341
#[2,] 0.04288691 0.007022297 -0.028842313
#[3,] 0.01058595 0.001696835 -0.007192284

#[[3]]
#            [,1]        [,2]         [,3]
#[1,]  0.27286729  0.09885234 -0.075162599
#[2,]  0.12991059  0.04437864 -0.041153317
#[3,] -0.01304611 -0.01009507 -0.007144035

library(RColorBrewer)
cs<-brewer.pal(n=11,"RdBu")
bs<-c(0.01,0.02,0.04,0.08,0.15,0.3)
brks<-c(-1 * rev(bs),bs)
library(fields)
pdf("ColorEpiPc1x2.pdf",width=7,height=7)
par(mfrow=c(2,2))
for(i in 1:3){
	image(gfit[[i]],breaks=brks,col=cs,axes=FALSE)
	title(main=paste("Fitness, PC1 = ",pcs[i],sep=""))
	axis(1,at=c(0,.5,1),c("aa","Aa","AA"))
	axis(2,at=c(0,.5,1),c("aa","Aa","AA"))
	box()
}
plot(c(0,1),c(0,1),type='n',axes=F,xlab="",ylab="")
image.plot(gfit[[i]],breaks=brks,col=cs,axes=FALSE,legend.only=TRUE,legend.lab="Rel. fitness")
dev.off()

library(RColorBrewer)
cl<-1.6
cm<-1.5
cn<-1.5
cs<-brewer.pal(n=3,"BrBG")
pdf("epiExample.pdf",width=6,height=12)
tits<-c("PC1 = -1","PC1 = 0","PC1 = 1")
par(mfrow=c(3,1))
par(mar=c(5,5,3,.5))
for(i in 1:3){
	barplot(gfit[[i]],beside=TRUE,xlab="Genotype",col=cs,names.arg=c("aa","Aa","AA"),ylab="Effect",cex.lab=cl,cex.names=cn)
	title(main=tits[i],cex.main=cm)
	if(i==2){legend(10,.05,c("aa","Aa","AA"),fill=cs,bty='n',cex=1.55)}
}
dev.off()

save(list=ls(),file="models.rdat")

## leave-one-out cross-validation
V<-15*2 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
	pcmod[,i]<-pcmod[,i] * pc_AC
}
dfAC<-data.frame(p_AC[,2],mod,pcmod)
Ni<-dim(dfAC)[1]
kgrp<-sample(1:Ni,Ni,replace=FALSE)
predAC<-rep(NA,Ni)
for(i in 101:200){
	dfACx<-dfAC[-c(kgrp==i),]
        o<-bms(dfACx,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
	j<-which(kgrp==i)
	predAC[j]<-sum(coef(o,order.by.pip=FALSE)[,2]*dfAC[j,-1])
}	
cor.test(dfAC[,1],predAC)
cov(dfAC[,1],predAC)/var(dfAC[,1])

V<-15*2 ## variables
form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
	pcmod[,i]<-pcmod[,i] * pc_MM
}
dfMM<-data.frame(p_MM[,2],mod,pcmod)
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

##### version without snp 6 but with PC1 for indel and for rest of mel stripe
pc<-prcomp(c_AC_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.9763 0.4129 0.37197 0.30890 0.26942
#Proportion of Variance 0.6665 0.1192 0.09676 0.06673 0.05076
#Cumulative Proportion  0.6665 0.7857 0.88251 0.94924 1.00000

# loadings
#             PC1        PC2         PC3          PC4        PC5
#[1,] -0.24204656  0.3938253 -0.83192866  0.193315555  0.2384092
#[2,]  0.02979171 -0.3573832 -0.05569658  0.885077280 -0.2914203
#[3,] -0.73888724  0.3315295  0.21356837 -0.007702461 -0.5463177
#[4,]  0.59561276  0.7023205  0.09707544  0.156509549 -0.3436168
#[5,]  0.19954200 -0.3376582 -0.49975744 -0.393335402 -0.6646049

pc_AC<-pc$x[,1]

pc<-prcomp(c_AC_ms,center=TRUE,scale=FALSE)
#Importance of components:
#                           PC1    PC2     PC3     PC4     PC5    PC6     PC7
#Standard deviation     1.24820 1.1286 1.06322 1.00518 0.95218 0.9070 0.81269
#Proportion of Variance 0.09394 0.0768 0.06816 0.06092 0.05467 0.0496 0.03982
#Cumulative Proportion  0.09394 0.1707 0.23891 0.29983 0.35450 0.4041 0.44392

ms_AC<-pc$x[,1]

cor(ms_AC,pc_AC)
#[1] -0.01612651

library(BMS)
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
bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
	bpipsAC[,i]<-oo[,1]
	bmavAC[,i]<-oo[,2]
}

paIds<-row.names(oo)

pc<-prcomp(c_MM_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.911 0.4188 0.36091 0.34237 0.24600
#Proportion of Variance 0.632 0.1335 0.09918 0.08925 0.04608
#Cumulative Proportion  0.632 0.7655 0.86467 0.95392 1.00000

# loadings
#             PC1        PC2         PC3        PC4         PC5
#[1,] -0.28343661  0.0572193 -0.90821797  0.2950037 -0.06710139
#[2,]  0.07555966 -0.4847130 -0.32289846 -0.8008655  0.11702633
#[3,] -0.75147812  0.3571962  0.13397790 -0.2733475  0.46370522
#[4,]  0.55736577  0.7061458 -0.22872329 -0.2411958  0.28321688
#[5,]  0.19641402 -0.3681689 -0.02474436  0.3724151  0.82859394

pc_MM<-pc$x[,1]

pc<-prcomp(c_MM_ms,center=TRUE,scale=FALSE)
#Importance of components:
#                           PC1    PC2     PC3     PC4     PC5    PC6     PC7
#Standard deviation     1.14247 1.05988 1.04498 1.02112 0.93840 0.91526 0.85776
#Proportion of Variance 0.08061 0.06937 0.06744 0.06439 0.05438 0.05173 0.04544
#Cumulative Proportion  0.08061 0.14998 0.21742 0.28181 0.33620 0.38793 0.43337

ms_MM<-pc$x[,1]

cor(ms_MM,pc_MM)
#[1] -0.1422148

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
bmsMM<-vector("list",5)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
	oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
	bpipsMM[,i]<-oo[,1]
	bmavMM[,i]<-oo[,2]
}

cs<-rep(c("cornsilk3","coral3","gold1","dodgerblue4","darkorange3","chartreuse4"),c(5,10,5,10,5,10))

pdf("ColorSurvPips5snpPcMs.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

pdf("ColorSurvModAvEffs5snpPcMs.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.2,.2),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Survival on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.2,.2),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Survival on MM",cex.main=1.5)
dev.off()

## number of additive effects, number of interactions in the model
et<-rep(c(1:6),c(5,10,5,10,5,10))
Nx<-10000
noAC<-matrix(NA,nrow=Nx,ncol=6)
noMM<-matrix(NA,nrow=Nx,ncol=6)
for(i in 1:Nx){
	sAC<-rbinom(n=45,size=1,prob=apply(bpipsAC,1,mean))
	sMM<-rbinom(n=45,size=1,prob=apply(bpipsMM,1,mean))
	noAC[i,]<-tapply(X=sAC,INDEX=et,sum)
	noMM[i,]<-tapply(X=sMM,INDEX=et,sum)
}

mnsd<-function(x=NA){
	c(mean(x),sd(x))
}

noQTL_AC<-matrix(NA,nrow=8,ncol=2)
noQTL_MM<-matrix(NA,nrow=8,ncol=2)
for(i in 1:6){
	noQTL_AC[i,]<-mnsd(noAC[,i])
	noQTL_MM[i,]<-mnsd(noMM[,i])
}
noQTL_AC[7,]<-mnsd(apply(noAC[,2:6],1,sum))
noQTL_MM[7,]<-mnsd(apply(noMM[,2:6],1,sum))
noQTL_AC[8,]<-mnsd(apply(noAC,1,sum))
noQTL_MM[8,]<-mnsd(apply(noMM,1,sum))

nms<-c("G","GxG","GxPC1","GxGxPC1","GxMs","GxGxMs","epis.","all")
pdf("ColorNoEffsMs.pdf",width=6,height=5)
par(mar=c(6,4.5,0.5,1))
cs<-c("orange","cadetblue")
ox<-barplot(t(as.matrix(cbind(noQTL_AC[,1],noQTL_MM[,1]))),ylim=c(0,8),col=cs,beside=TRUE,xlab="",ylab="No. of effects",cex.lab=1.5,las=2,cex.names=1.2)
noQTL<-cbind(noQTL_AC,noQTL_MM)
oxm<-apply(ox,2,mean)
axis(1,at=oxm,nms,las=2,cex.axis=1.2)
ox<-as.vector(ox)
ox<-ox[c(1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16)]
segments(ox,noQTL[,1]+noQTL[,2],ox,noQTL[,1]-noQTL[,2])
legend(1,5.8,c("AC","MM"),fill=cs,bty='n',cex=1.2)
box()
dev.off()


##############################################################################
## fits of 'top' model and models with pips > 0.5... these are the same
xx<-which(apply(bpipsAC,1,mean) > .5)
## 6 8
o<-glm(dfAC[,1] ~ as.matrix(dfAC[,xx+1]),family="binomial")
summary(o)
#Coefficients:
#                               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                     -1.4391     0.2011  -7.157 8.24e-13 ***
#as.matrix(dfAC[, xx + 1])V6      2.1053     0.7271   2.896  0.00378 ** 
#as.matrix(dfAC[, xx + 1])V1.V3  -1.5453     0.5209  -2.967  0.00301 ** 
#    Null deviance: 240.71  on 215  degrees of freedom
#Residual deviance: 224.95  on 213  degrees of freedom
#AIC: 230.95

## Efron's pseudo r2
mn<-mean(dfAC[,1])
1-(sum((dfAC[,1]-o$fitted.values)^2)/sum((dfAC[,1]-mn)^2))
#[1] 0.08097351

## comparing to case where all terms with pp > 0.1 are included
xx<-which(apply(bpipsAC,1,mean) > .1)
##  4  6  7  8  9 10 13 14 15 18 21
o<-glm(dfAC[,1] ~ as.matrix(dfAC[,xx+1]),family="binomial")
summary(o)
#Coefficients:
#                               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                     -1.5309     0.2917  -5.249 1.53e-07 ***
#as.matrix(dfAC[, xx + 1])V4     -0.8559     0.3808  -2.247   0.0246 *  
#as.matrix(dfAC[, xx + 1])V6      5.2336     3.4741   1.506   0.1320    
#as.matrix(dfAC[, xx + 1])V1.V2   1.4741     1.6857   0.874   0.3819    
#as.matrix(dfAC[, xx + 1])V1.V3  -2.0473     0.9999  -2.047   0.0406 *  
#as.matrix(dfAC[, xx + 1])V1.V4   0.1972     1.1848   0.166   0.8678    
#as.matrix(dfAC[, xx + 1])V1.V5  -1.7193     1.1768  -1.461   0.1440    
#as.matrix(dfAC[, xx + 1])V2.V4   0.6583     1.0109   0.651   0.5149    
#as.matrix(dfAC[, xx + 1])V2.V5  -6.0606     2.8971  -2.092   0.0364 *  
#as.matrix(dfAC[, xx + 1])V2.V6  -3.3460     1.7971  -1.862   0.0626 .  
#as.matrix(dfAC[, xx + 1])V3.V6  11.2613     7.2973   1.543   0.1228    
#as.matrix(dfAC[, xx + 1])V5.V6  -0.2194    13.3376  -0.016   0.9869    
#    Null deviance: 240.71  on 215  degrees of freedom
#Residual deviance: 207.22  on 204  degrees of freedom
#AIC: 231.22
mn<-mean(dfAC[,1])
1-(sum((dfAC[,1]-o$fitted.values)^2)/sum((dfAC[,1]-mn)^2))
## [1] 0.1569219
########################################################################################################################

fit_AC<-read.table("exfit_AC.txt",header=TRUE)
fit_MM<-read.table("exfit_MM.txt",header=TRUE)

##### version without snp 6 but with PC1
pc<-prcomp(c_AC_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.9763 0.4129 0.37197 0.30890 0.26942
#Proportion of Variance 0.6665 0.1192 0.09676 0.06673 0.05076
#Cumulative Proportion  0.6665 0.7857 0.88251 0.94924 1.00000

pc_AC<-pc$x[,1]

library(BMS)
V<-15*2 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
        pcmod[,i]<-pcmod[,i] * pc_AC
}
dfAC<-data.frame(fit_AC$exfit,blk2=as.numeric(fit_AC$blk==2),blk3=as.numeric(fit_AC$blk==3),mod,pcmod)
#lexfit<-log(fit_AC$exfit/(1-fit_AC$exfit))
#dfAC<-data.frame(lexfit,mod,pcmod)

bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
        bpipsAC[,i]<-oo[-c(1:2),1]
        bmavAC[,i]<-oo[-c(1:2),2]
}

paIds<-row.names(oo)[-c(1:2)]

pc<-prcomp(c_MM_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.911 0.4188 0.36091 0.34237 0.24600
#Proportion of Variance 0.632 0.1335 0.09918 0.08925 0.04608
#Cumulative Proportion  0.632 0.7655 0.86467 0.95392 1.00000

pc_MM<-pc$x[,1]

form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
        pcmod[,i]<-pcmod[,i] * pc_MM
}
dfMM<-data.frame(fit_MM$exfit,blk2=as.numeric(fit_MM$blk==2),blk3=as.numeric(fit_MM$blk==3),mod,pcmod)
#lexfit<-log(fit_MM$exfit/(1-fit_MM$exfit))
#dfMM<-data.frame(lexfit,mod,pcmod)
bmsMM<-vector("list",5)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
        bpipsMM[,i]<-oo[-c(1:2),1]
        bmavMM[,i]<-oo[-c(1:2),2]
}

cs<-rep(c("cornsilk3","coral3","gold1","dodgerblue4"),c(5,10,5,10))

pdf("ColorExFitPips5snpPC.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

pdf("ColorExFitModAvEffs5snpPC.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

## number of additive effects, number of interactions in the model
et<-rep(c(1:4),c(5,10,5,10))
Nx<-10000
noAC<-matrix(NA,nrow=Nx,ncol=4)
noMM<-matrix(NA,nrow=Nx,ncol=4)
for(i in 1:Nx){
        sAC<-rbinom(n=30,size=1,prob=apply(bpipsAC,1,mean))
        sMM<-rbinom(n=30,size=1,prob=apply(bpipsMM,1,mean))
        noAC[i,]<-tapply(X=sAC,INDEX=et,sum)
        noMM[i,]<-tapply(X=sMM,INDEX=et,sum)
}

mnsd<-function(x=NA){
        c(mean(x),sd(x))
}

noQTL_AC<-matrix(NA,nrow=6,ncol=2)
noQTL_MM<-matrix(NA,nrow=6,ncol=2)
for(i in 1:4){
        noQTL_AC[i,]<-mnsd(noAC[,i])
        noQTL_MM[i,]<-mnsd(noMM[,i])
}
noQTL_AC[5,]<-mnsd(apply(noAC[,2:4],1,sum))
noQTL_MM[5,]<-mnsd(apply(noMM[,2:4],1,sum))
noQTL_AC[6,]<-mnsd(apply(noAC,1,sum))
noQTL_MM[6,]<-mnsd(apply(noMM,1,sum))

nms<-c("G","GxG","GxPC1","GxGxPC1","epis.","all")
pdf("ColorNoEffsExFit.pdf",width=6,height=5)
par(mar=c(6,4.5,0.5,1))
cs<-c("orange","cadetblue")
ox<-barplot(t(as.matrix(cbind(noQTL_AC[,1],noQTL_MM[,1]))),ylim=c(0,11),col=cs,beside=TRUE,xlab="",ylab="No. of effects",cex.lab=1.5,las=2,cex.names=1.2)
noQTL<-rbind(noQTL_AC,noQTL_MM)
oxm<-apply(ox,2,mean)
axis(1,at=oxm,nms,las=2,cex.axis=1.2)
ox<-as.vector(ox)
ox<-ox[c(1,3,5,7,9,11,2,4,6,8,10,12)]
segments(ox,noQTL[,1]+noQTL[,2],ox,noQTL[,1]-noQTL[,2])
legend(1,10,c("AC","MM"),fill=cs,bty='n',cex=1.2)
box()
dev.off()


######## 3-way epistasis ##############
fit_AC<-read.table("exfit_AC.txt",header=TRUE)
fit_MM<-read.table("exfit_MM.txt",header=TRUE)

##### version without snp 6 but with PC1
pc<-prcomp(c_AC_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.9763 0.4129 0.37197 0.30890 0.26942
#Proportion of Variance 0.6665 0.1192 0.09676 0.06673 0.05076
#Cumulative Proportion  0.6665 0.7857 0.88251 0.94924 1.00000

pcAC_ld<-pc$rotation
#             PC1        PC2         PC3          PC4        PC5
#[1,] -0.24204656  0.3938253 -0.83192866  0.193315555  0.2384092
#[2,]  0.02979171 -0.3573832 -0.05569658  0.885077280 -0.2914203
#[3,] -0.73888724  0.3315295  0.21356837 -0.007702461 -0.5463177
#[4,]  0.59561276  0.7023205  0.09707544  0.156509549 -0.3436168
#[5,]  0.19954200 -0.3376582 -0.49975744 -0.393335402 -0.6646049


library(BMS)
V<-19 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
by3<-mod[,1:4]
by3[,1]<-c_AC_pip[,3] * c_AC_pip[,1] * c_AC_pip[,4]
by3[,2]<-c_AC_pip[,3] * c_AC_pip[,1] * c_AC_pip[,5]
by3[,3]<-c_AC_pip[,3] * c_AC_pip[,4] * c_AC_pip[,5]
by3[,4]<-c_AC_pip[,1] * c_AC_pip[,4] * c_AC_pip[,5]

dfAC<-data.frame(fit_AC$exfit,blk2=as.numeric(fit_AC$blk==2),blk3=as.numeric(fit_AC$blk==3),mod,by3)
#lexfit<-log(fit_AC$exfit/(1-fit_AC$exfit))
#dfAC<-data.frame(lexfit,mod,pcmod)

bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
        bpipsAC[,i]<-oo[-c(1:2),1]
        bmavAC[,i]<-oo[-c(1:2),2]
}

V<-19 ## variables
form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
by3<-mod[,1:4]
by3[,1]<-c_MM_pip[,3] * c_MM_pip[,1] * c_MM_pip[,4]
by3[,2]<-c_MM_pip[,3] * c_MM_pip[,1] * c_MM_pip[,5]
by3[,3]<-c_MM_pip[,3] * c_MM_pip[,4] * c_MM_pip[,5]
by3[,4]<-c_MM_pip[,1] * c_MM_pip[,4] * c_MM_pip[,5]

dfMM<-data.frame(fit_MM$exfit,blk2=as.numeric(fit_MM$blk==2),blk3=as.numeric(fit_MM$blk==3),mod,by3)
#lexfit<-log(fit_AC$exfit/(1-fit_AC$exfit))
#dfAC<-data.frame(lexfit,mod,pcmod)

bmsMM<-vector("list",5)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
        bpipsMM[,i]<-oo[-c(1:2),1]
        bmavMM[,i]<-oo[-c(1:2),2]
}


cs<-rep(c("cornsilk3","coral3","gold1"),c(5,10,4))
paIds<-row.names(oo)[-c(1:2)]

paIds[16:19]<-c("V3.V1.V4","V3.V1.V5","V3.V4.V5","V1.V4.V5")


pdf("ColorExFitPips5snp3way.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

pdf("ColorExFitModAvEffs5snp3way.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

######## all the epistasis ##############
fit_AC<-read.table("exfit_AC.txt",header=TRUE)
fit_MM<-read.table("exfit_MM.txt",header=TRUE)


library(BMS)
V<-31 ## variables
form<-p_AC[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])

dfAC<-data.frame(fit_AC$exfit,blk2=as.numeric(fit_AC$blk==2),blk3=as.numeric(fit_AC$blk==3),mod)

bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
        bpipsAC[,i]<-oo[-c(1:2),1]
        bmavAC[,i]<-oo[-c(1:2),2]
}

V<-31 ## variables
form<-p_MM[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])

dfMM<-data.frame(fit_MM$exfit,blk2=as.numeric(fit_MM$blk==2),blk3=as.numeric(fit_MM$blk==3),mod)

bmsMM<-vector("list",5)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
        bpipsMM[,i]<-oo[-c(1:2),1]
        bmavMM[,i]<-oo[-c(1:2),2]
}

paIds<-row.names(oo)[-c(1:2)]

cs<-rep(c("cornsilk3","coral3","gold1","dodgerblue4","black"),c(5,10,10,5,1))
paIds<-row.names(oo)[-c(1:2)]

pdf("ColorExFitPips5snpAllEpi.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

pdf("ColorExFitModAvEffs5snpAllEpi.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

## number of additive, 2-way, 3-way, 4-way and 5-way 
et<-rep(c(1:5),c(5,10,10,5,1))
Nx<-10000
noAC<-matrix(NA,nrow=Nx,ncol=5)
noMM<-matrix(NA,nrow=Nx,ncol=5)
for(i in 1:Nx){
        sAC<-rbinom(n=31,size=1,prob=apply(bpipsAC,1,mean))
        sMM<-rbinom(n=31,size=1,prob=apply(bpipsMM,1,mean))
        noAC[i,]<-tapply(X=sAC,INDEX=et,sum)
        noMM[i,]<-tapply(X=sMM,INDEX=et,sum)
}

mnsd<-function(x=NA){
        c(mean(x),sd(x))
}

noQTL_AC<-matrix(NA,nrow=7,ncol=2)
noQTL_MM<-matrix(NA,nrow=7,ncol=2)
for(i in 1:5){
        noQTL_AC[i,]<-mnsd(noAC[,i])
        noQTL_MM[i,]<-mnsd(noMM[,i])
}
noQTL_AC[6,]<-mnsd(apply(noAC[,2:5],1,sum))
noQTL_MM[6,]<-mnsd(apply(noMM[,2:5],1,sum))
noQTL_AC[7,]<-mnsd(apply(noAC,1,sum))
noQTL_MM[7,]<-mnsd(apply(noMM,1,sum))

nms<-c("add.","2-way","3-way","4-way","5-way","epis.","all")
pdf("ColorNoEffsExFitEpi.pdf",width=6,height=5)
par(mar=c(6,4.5,0.5,1))
cs<-c("orange","cadetblue")
ox<-barplot(t(as.matrix(cbind(noQTL_AC[,1],noQTL_MM[,1]))),ylim=c(0,10.5),col=cs,beside=TRUE,xlab="",ylab="No. of effects",cex.lab=1.5,las=2,cex.names=1.2)
noQTL<-rbind(noQTL_AC,noQTL_MM)
oxm<-apply(ox,2,mean)
axis(1,at=oxm,nms,las=2,cex.axis=1.2)
ox<-as.vector(ox)
ox<-ox[c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]
segments(ox,noQTL[,1]+noQTL[,2],ox,noQTL[,1]-noQTL[,2])
legend(1,10,c("AC","MM"),fill=cs,bty='n',cex=1.2)
box()
dev.off()

## pairwise LD for pip SNPs
LD<-(cor(rbind(c_AC_pip[,-6],c_MM_pip[,-6])))^2
pdf("indelPipLD.pdf",width=5,height=5)
par(mar=c(4.5,4.5,0.5,.5))
plot(sort(LD[upper.tri(LD)]),pch=19,cex.lab=1.4,cex.axis=1.1,xlab="SNP pair",ylab="LD (r2)")
mtext(paste("mean = ",round(mean(LD[upper.tri(LD)]),3),sep=""),3,adj=0.1,line=-2)
dev.off()

## bayesian ridge
fit_AC<-read.table("exfit_AC.txt",header=TRUE)
fit_MM<-read.table("exfit_MM.txt",header=TRUE)

library(rstan)

modelStringRidge<-"data{
    int N_train;             // # training observations
    int N_test;              // # test observations
    int N_pred;              // # predictor variables
    vector[N_train] y_train; // training outcomes
    matrix[N_train, N_pred] X_train; // training data
    matrix[N_test, N_pred] X_test;   // testing data
}
parameters{
    real<lower=0> sigma;   // error SD
    real<lower=0> sigma_B; // hierarchical SD across betas
    vector[N_pred] beta;   // regression beta weights
}
model{
  // group-level (hierarchical) SD across betas
  sigma_B ~ cauchy(0, 1);

  // model error SD
  sigma ~ normal(0, 1);

  // beta prior (provides ridge regularization)
  beta ~ normal(0, sigma_B);

  // model likelihood
    y_train ~ normal(X_train*beta, sigma);
}
generated quantities{
    real y_test[N_test]; // test data predictions
    for(i in 1:N_test){
        y_test[i] = normal_rng(X_test[i,] * beta, sigma);
    }
}
"

stanDsoRidge<-stan_model(model_code=modelStringRidge)

modelStringLasso<-"data{
    int N_train;             // # training observations
    int N_test;              // # test observations
    int N_pred;              // # predictor variables
    vector[N_train] y_train; // training outcomes
    matrix[N_train, N_pred] X_train; // training data
    matrix[N_test, N_pred] X_test;   // testing data
}
parameters{
    real<lower=0> sigma;   // error SD
    real<lower=0> sigma_B; // hierarchical SD across betas
    vector[N_pred] beta;   // regression beta weights
}
model{
  // group-level (hierarchical) SD across betas
  sigma_B ~ cauchy(0, 1);

  // model error SD
  sigma ~ normal(0, 1);

  // beta prior (provides lasso regularization)
  beta ~ double_exponential(0, sigma_B);

  // model likelihood
    y_train ~ normal(X_train*beta, sigma);
}
generated quantities{
    real y_test[N_test]; // test data predictions
    for(i in 1:N_test){
        y_test[i] = normal_rng(X_test[i,] * beta, sigma);
    }
}
"

stanDsoLasso<-stan_model(model_code=modelStringLasso)

V<-31 ## variables
form<-p_AC[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])

X<-data.frame(as.numeric(fit_AC$blk==1),as.numeric(fit_AC$blk==2),as.numeric(fit_AC$blk==3),
	mod)

nas<-which(is.na(fit_AC$exfit)==TRUE)
kp<-1:200
kp<-kp[-nas]

sl_AC<-list(N_train=197, N_test=19, N_pred=V+3,
	y_train=fit_AC$exfit[kp], X_train=X[kp,],X_test=X[201:219,])

fit_ridge_AC <- sampling(stanDsoRidge, sl_AC, iter = 2000, 
                            warmup = 500, chains = 3, cores = 3)
fit_lasso_AC <- sampling(stanDsoLasso, sl_AC, iter = 2000, 
                            warmup = 500, chains = 3, cores = 3)


bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
        bpipsAC[,i]<-oo[-c(1:2),1]
        bmavAC[,i]<-oo[-c(1:2),2]
}

V<-31 ## variables
form<-p_MM[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])

dfMM<-data.frame(fit_MM$exfit,blk2=as.numeric(fit_MM$blk==2),blk3=as.numeric(fit_MM$blk==3),mod)

bmsMM<-vector("list",5)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

######## all the epistasis and dominance##############
fit_AC<-read.table("exfit_AC.txt",header=TRUE)
fit_MM<-read.table("exfit_MM.txt",header=TRUE)


library(BMS)
V<-33 ## variables
form<-p_AC[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
d3<-1-abs(g_AC_pip[,3]-1)## het = 1, both homos = 0
d4<-1-abs(g_AC_pip[,4]-1)

dfAC<-data.frame(fit_AC$exfit,blk2=as.numeric(fit_AC$blk==2),blk3=as.numeric(fit_AC$blk==3),mod,d3,d4)

bmsAC<-vector("list",5)
for(i in 1:5){
        bmsAC[[i]]<-bms(dfAC,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsAC<-matrix(NA,nrow=V,ncol=5)
bmavAC<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsAC[[i]],order.by.pip=FALSE)
        bpipsAC[,i]<-oo[-c(1:2),1]
        bmavAC[,i]<-oo[-c(1:2),2]
}

V<-33 ## variables
form<-p_MM[,2] ~ .^5
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
d3<-1-abs(g_MM_pip[,3]-1)## het = 1, both homos = 0
d4<-1-abs(g_MM_pip[,4]-1)

dfMM<-data.frame(fit_MM$exfit,blk2=as.numeric(fit_MM$blk==2),blk3=as.numeric(fit_MM$blk==3),mod,d3,d4)

bmsMM<-vector("list",5)
for(i in 1:5){
        bmsMM[[i]]<-bms(dfMM,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE,fixed.reg=c(1,2))
}

bpipsMM<-matrix(NA,nrow=V,ncol=5)
bmavMM<-matrix(NA,nrow=V,ncol=5)
for(i in 1:5){
        oo<-coef(bmsMM[[i]],order.by.pip=FALSE)
        bpipsMM[,i]<-oo[-c(1:2),1]
        bmavMM[,i]<-oo[-c(1:2),2]
}

paIds<-row.names(oo)[-c(1:2)]

cs<-rep(c("cornsilk3","coral3","gold1","dodgerblue4","black","gray"),c(5,10,10,5,1,2))
paIds<-row.names(oo)[-c(1:2)]

pdf("ColorExFitPips5snpAllEpiDom.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bpipsAC,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bpipsMM,1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

pdf("ColorExFitModAvEffs5snpAllEpiDom.pdf",width=7,height=9)
par(mfrow=c(2,1))
par(mar=c(5.75,4.5,3,0.5))
barplot(apply(bmavAC,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(A) Expected fitness on A/C",cex.main=1.5)
barplot(apply(bmavMM,1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
title(main="(B) Expected fitness on MM",cex.main=1.5)
dev.off()

## number of additive, 2-way, 3-way, 4-way and 5-way 
et<-rep(c(1:6),c(5,10,10,5,1,2))
Nx<-10000
noAC<-matrix(NA,nrow=Nx,ncol=6)
noMM<-matrix(NA,nrow=Nx,ncol=6)
for(i in 1:Nx){
        sAC<-rbinom(n=33,size=1,prob=apply(bpipsAC,1,mean))
        sMM<-rbinom(n=33,size=1,prob=apply(bpipsMM,1,mean))
        noAC[i,]<-tapply(X=sAC,INDEX=et,sum)
        noMM[i,]<-tapply(X=sMM,INDEX=et,sum)
}

mnsd<-function(x=NA){
        c(mean(x),sd(x))
}

noQTL_AC<-matrix(NA,nrow=8,ncol=2)
noQTL_MM<-matrix(NA,nrow=8,ncol=2)
for(i in 1:6){
        noQTL_AC[i,]<-mnsd(noAC[,i])
        noQTL_MM[i,]<-mnsd(noMM[,i])
}
noQTL_AC[7,]<-mnsd(apply(noAC[,2:5],1,sum))
noQTL_MM[7,]<-mnsd(apply(noMM[,2:5],1,sum))
noQTL_AC[8,]<-mnsd(apply(noAC,1,sum))
noQTL_MM[8,]<-mnsd(apply(noMM,1,sum))

nms<-c("add.","2-way","3-way","4-way","5-way","dom","epis.","all")
pdf("ColorNoEffsExFitEpiDom.pdf",width=6,height=5)
par(mar=c(6,4.5,0.5,1))
cs<-c("orange","cadetblue")
ox<-barplot(t(as.matrix(cbind(noQTL_AC[,1],noQTL_MM[,1]))),ylim=c(0,10.5),col=cs,beside=TRUE,xlab="",ylab="No. of effects",cex.lab=1.5,las=2,cex.names=1.2)
noQTL<-rbind(noQTL_AC,noQTL_MM)
oxm<-apply(ox,2,mean)
axis(1,at=oxm,nms,las=2,cex.axis=1.2)
ox<-as.vector(ox)
ox<-ox[c(1,3,5,7,9,11,13,15,2,4,6,8,10,12,14,16)]
segments(ox,noQTL[,1]+noQTL[,2],ox,noQTL[,1]-noQTL[,2])
legend(1,10,c("AC","MM"),fill=cs,bty='n',cex=1.2)
box()
dev.off()



#############################################################################################
## fits by block

##### version without snp 6 but with PC1
pc<-prcomp(c_AC_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.9763 0.4129 0.37197 0.30890 0.26942
#Proportion of Variance 0.6665 0.1192 0.09676 0.06673 0.05076
#Cumulative Proportion  0.6665 0.7857 0.88251 0.94924 1.00000

pc_AC<-pc$x[,1]

library(BMS)
V<-15*2 ## variables
form<-p_AC[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_AC_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
        pcmod[,i]<-pcmod[,i] * pc_AC
}

dfAC<-data.frame(fit_AC$exfit,mod,pcmod)
bpipsACblk<-vector("list",3)
bmavACblk<-vector("list",3)
for(k in 1:3){
	kblk<-which(fit_AC$blk==k)
	dfACblk<-dfAC[kblk,]

	bmsACblk<-vector("list",5)
	for(i in 1:5){
	       	bmsACblk[[i]]<-bms(dfACblk,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
	}

	bpipsACblk[[k]]<-matrix(NA,nrow=V,ncol=5)
	bmavACblk[[k]]<-matrix(NA,nrow=V,ncol=5)
	for(i in 1:5){
	        oo<-coef(bmsACblk[[i]],order.by.pip=FALSE)
       		bpipsACblk[[k]][,i]<-oo[,1]
        	bmavACblk[[k]][,i]<-oo[,2]
	}
}

paIds<-row.names(oo)[-c(1:2)]

pc<-prcomp(c_MM_pip[,-6],center=TRUE,scale=FALSE)
#Importance of components:
#                          PC1    PC2     PC3     PC4     PC5
#Standard deviation     0.911 0.4188 0.36091 0.34237 0.24600
#Proportion of Variance 0.632 0.1335 0.09918 0.08925 0.04608
#Cumulative Proportion  0.632 0.7655 0.86467 0.95392 1.00000

pc_MM<-pc$x[,1]

form<-p_MM[,2] ~ .^2
mod<-model.matrix(form,dat=as.data.frame(c_MM_pip[,-6]))
mod<-cbind(mod[,-1])
pcmod<-mod
for(i in 1:(V/2)){
        pcmod[,i]<-pcmod[,i] * pc_MM
}

dfMM<-data.frame(fit_MM$exfit,mod,pcmod)
bpipsMMblk<-vector("list",3)
bmavMMblk<-vector("list",3)
for(k in 1:3){
	kblk<-which(fit_MM$blk==k)
	dfMMblk<-dfMM[kblk,]

	bmsMMblk<-vector("list",5)
	for(i in 1:5){
	       	bmsMMblk[[i]]<-bms(dfMMblk,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
	}

	bpipsMMblk[[k]]<-matrix(NA,nrow=V,ncol=5)
	bmavMMblk[[k]]<-matrix(NA,nrow=V,ncol=5)
	for(i in 1:5){
	        oo<-coef(bmsMMblk[[i]],order.by.pip=FALSE)
       		bpipsMMblk[[k]][,i]<-oo[,1]
        	bmavMMblk[[k]][,i]<-oo[,2]
	}
}

cs<-rep(c("cornsilk3","coral3","gold1","dodgerblue4"),c(5,10,5,10))

pdf("ColorExFitPips5snpPC_blk.pdf",width=12,height=12)
par(mfrow=c(3,2))
par(mar=c(5.75,4.5,3,0.5))
for(k in 1:3){
	barplot(apply(bpipsACblk[[k]],1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
	title(main=paste("(",LETTERS[k],") Expected fitness on A/C, blk =",k,sep=""),cex.main=1.5)
	barplot(apply(bpipsMMblk[[k]],1,mean),ylim=c(0,1),names.arg=paIds,col=cs,xlab="Coefficient",ylab="PIP",cex.lab=1.4,las=2,cex.names=.7)
	title(main=paste("(",LETTERS[3+k],") Expected fitness on MM, blk =",k,sep=""),cex.main=1.5)
}
dev.off()

pdf("ColorExFitModAvEffs5snpPC_blk.pdf",width=12,height=12)
par(mfrow=c(3,2))
par(mar=c(5.75,4.5,3,0.5))
for(k in 1:3){
	barplot(apply(bmavACblk[[k]],1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
	title(main=paste("(",LETTERS[k],") Expected fitness on A/C, blk =",k,sep=""),cex.main=1.5)
	barplot(apply(bmavMMblk[[k]],1,mean),ylim=c(-.5,.5),names.arg=paIds,col=cs,xlab="Coefficient",ylab="Effect",cex.lab=1.4,las=2,cex.names=.7)
	title(main=paste("(",LETTERS[3+k],") Expected fitness on MM, blk =",k,sep=""),cex.main=1.5)
}
dev.off()

## not sure I trust the block ones, at least not yet!
