### Load in the R libraries ###
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)

### Load in functions to make QQ-plot plots ###
source("MAPIT/QQPlot.R")

#NOTE: This code assumes that the basic C++ functions are set up on the computer in use. If not, the MAPIT functions and Rcpp packages will not work properly. Mac users please refer to the homebrew applications and install the gcc commands listed in the README.md file before running the rest of the code [Warning: This step may take about an hour...].

### Load in the C++ MAPIT wrapper functions ###
source("MAPIT/Standard Version/MAPIT.R"); sourceCpp("MAPIT/Standard Version/MAPIT.cpp")

## read in color data
cdat<-read.table("../Gemma_mel/pheno_color.txt",header=FALSE)
## 437 inds x 2 traits

## read in genotype data
gdat<-as.matrix(read.table("../Gemma_mel/g_tchum.txt",header=FALSE))

## subset mel stripe
snps<-as.matrix(read.table("../Gemma_mel/tchumSnpTable.txt",header=FALSE))
ms<-c(which(snps[,2]==128 & snps[,3] <6414835),which(snps[,2]==702.1 & snps[,3] < 4139489))
msSnps<-snps[ms,]

## standardize genotypes
X<-gdat[ms,]
Xmn<-apply(X,1,mean,na.rm=TRUE)
Xsd<-apply(X,1,sd,na.rm=TRUE)
L<-length(ms)
for(i in 1:L){
	X[i,]<-(X[i,]-Xmn[i])/Xsd[i]
}
msSnps<-msSnps[-c(which(is.na(apply(X,1,mean))==T)),]
drpSnps<-ms[c(which(is.na(apply(X,1,mean))==T))]

X<-X[-c(which(is.na(apply(X,1,mean))==T)),]## drop 2 invariant snps

yrg<-(cdat[,1]-mean(cdat[,1],na.rm=TRUE))/sd(cdat[,1],na.rm=TRUE)
ygb<-(cdat[,2]-mean(cdat[,2],na.rm=TRUE))/sd(cdat[,2],na.rm=TRUE)

mis<-which(is.na(yrg)==TRUE)

yrg<-yrg[-mis]
ygb<-ygb[-mis]
X<-X[,-mis]

## run mapit
mapit_rg<-MAPIT(X,as.matrix(yrg))
mapit_gb<-MAPIT(X,as.matrix(ygb))

#pdf("mapitEpi.pdf",width=5,height=5)
#par(mar=c(4,4,.7,.7))
#plot(-1 *log10(mapit_rg$pvalues),-1 *log10(mapit_gb$pvalues),xlab="-log10(p) RG",ylab="-log10(p) GB",
#	cex.lab=1.4,pch=19)
#gws<- -1 * log10(0.05/158)
#abline(v=gws,lty=2)
#abline(h=gws,lty=2)
#dev.off()

## PVE for sig SNPs
sum(mapit_rg$pves[mapit_rg$pvalues < (0.05/158)],na.rm=TRUE)
#[1] 3.100179
sum(mapit_gb$pves[mapit_gb$pvalues < (0.05/158)],na.rm=TRUE)
#[1] 0.5003422

## grab epistasis terms
gws<- -1 * log10(0.05/158)

sig<-which(mapit_rg$pvalues < (0.05/158) | mapit_gb$pvalues < (0.05/158))
eSnps<-msSnps[sig,]
#     V1  V2      V3
#[1,]  8 128 5359496
#[2,]  8 128 5690726
#[3,]  8 128 5690741
#[4,]  8 128 5746977
#[5,]  8 128 5746990


Gepi<-matrix(NA,nrow=5,ncol=437)
for(i in 1:5){
	a<-which(snps[,2] == eSnps[i,2] & snps[,3] == eSnps[i,3])
	Gepi[i,]<-gdat[a,]
	Gepi[i,]<-Gepi[i,]-mean(Gepi[i,])
}
eGmat<-matrix(NA,nrow=10,ncol=437)
k<-1
for(i in 1:4){for(j in (i+1):5){
	eGmat[k,]<-Gepi[i,] * Gepi[j,]
	k<-k+1
}}

Gout<-as.matrix(rbind(gdat,eGmat))
write.table(Gout,file="../Gemma_mel/g_epi_tchum.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

