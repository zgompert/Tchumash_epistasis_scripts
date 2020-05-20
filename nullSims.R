# fit BMS models for random sets of SNPs

library(BMS)
load("models.rdat")

nsims<-100

## sample five SNPs with similar allele freqs. to the color-associated SNPs
samsnps<-function(gen=g_AC,p=NA){
	px<-apply(gen,1,mean)/2
	sam<-rep(NA,5)
	# convert to maf
	px[px > 0.5]<-1-px[px > 0.5]
	p[p > 0.5]<-1-p[p > 0.5]
	good<-FALSE
	while(good==FALSE){
		ps<-rep(NA,5)
		for(i in 1:5){
			ps<-which(px >= (p[i]-0.025) & px <= (p[i]+0.025))
			sam[i]<-sample(ps,1,replace=FALSE)
		}
		if(length(unique(sam))==5){
			good<-TRUE
		}
	}
	gs<-t(gen[sam,])
	return(gs)	
}

## fit bms for one set of sampled snps
onesim<-function(fit=NA,c_gen=NA){
	c_gen<-as.matrix(c_gen)
	NN<-dim(c_gen)[[1]]
	V<-31 ## variables
	form<-rep(1,NN) ~ .^5
	mod<-model.matrix(form,dat=as.data.frame(c_gen))
	mod<-cbind(mod[,-1])
	df<-data.frame(fit$exfit,blk2=as.numeric(fit$blk==2),blk3=as.numeric(fit$blk==3),mod)
	obms<-bms(df,burn=10000,iter=200000,mprior="uniform",g="UIP",user.int=FALSE)
	oo<-coef(obms,order.by.pip=FALSE)
	bpips<-oo[-c(1,2),1]
	bmav<-oo[-c(1,2),2]
	## estiamte no. QTL
	et<-rep(c(1:5),c(5,10,10,5,1))
	Nx<-1000
	no<-matrix(NA,nrow=Nx,ncol=5)
	for(i in 1:Nx){
        	sa<-rbinom(n=31,size=1,prob=bpips)
		no[i,]<-tapply(X=sa,INDEX=et,sum)
	}
	noQTL<-rep(NA,7)
	mnEff<-rep(NA,7)
	for(i in 1:5){
		noQTL[i]<-mean(no[,i])
		mnEff[i]<-mean(abs(bmav[et==i]),na.rm=TRUE)
	}
	noQTL[6]<-mean(apply(no[,2:5],1,sum))
	noQTL[7]<-mean(apply(no,1,sum))
	mnEff[6]<-mean(abs(bmav[et > 1]),na.rm=TRUE)
	mnEff[7]<-mean(abs(bmav))
	c(noQTL,mnEff)
}

## get allele freqs.
ppac<-apply(g_AC_pip,2,mean)/2
ppmm<-apply(g_MM_pip,2,mean)/2

ppac<-ppac[1:5]
ppmm<-ppmm[1:5]

resAC<-matrix(NA,nrow=nsims,ncol=14)
resMM<-matrix(NA,nrow=nsims,ncol=14)

for(k in 1:nsims){
	cat(k,"\n")
	## AC
	## sample snps, igrnoe non LG snps
	g_sam<-samsnps(g_AC[snps[,1] > 0,],p=ppac)
	## center snps
	c_sam<-g_sam
	for(i in 1:5){
		c_sam[,i]<-c_sam[,i]-mean(c_sam[,i])
	}
	## fit model
	resAC[k,]<-onesim(fit=fit_AC,c_gen=c_sam)
	
	## MM
	## sample snps, igrnoe non LG snps
	g_sam<-samsnps(g_MM[snps[,1] > 0,],p=ppmm)
	## center snps
	c_sam<-g_sam
	for(i in 1:5){
		c_sam[,i]<-c_sam[,i]-mean(c_sam[,i])
	}
	## fit model
	resMM[k,]<-onesim(fit=fit_MM,c_gen=c_sam)
}

save(list=ls(),file="null.rdat")
write.table(file="nullAC.txt",resAC,row.names=FALSE,col.names=FALSE)
write.table(file="nullMM.txt",resMM,row.names=FALSE,col.names=FALSE)

noQTL_nullcomp_AC<-matrix(NA,nrow=7,ncol=2)
for(i in 1:7){
	noQTL_nullcomp_AC[i,1]<-noQTL_AC[i,1]/mean(resAC[,i])
	noQTL_nullcomp_AC[i,2]<-mean(resAC[,i] >= noQTL_AC[i,1])
}

noQTL_nullcomp_MM<-matrix(NA,nrow=7,ncol=2)
for(i in 1:7){
	noQTL_nullcomp_MM[i,1]<-noQTL_MM[i,1]/mean(resMM[,i])
	noQTL_nullcomp_MM[i,2]<-mean(resMM[,i] >= noQTL_MM[i,1])
}
## observed, x-fold and p-value
round(cbind(noQTL_MM[,1],noQTL_nullcomp_MM),2)
#     [,1] [,2] [,3]
#[1,] 0.84 1.25 0.20
#[2,] 1.46 1.14 0.24
#[3,] 1.49 1.20 0.22
#[4,] 0.83 1.37 0.19
#[5,] 0.17 1.62 0.11
#[6,] 3.96 1.22 0.21
#[7,] 4.81 1.23 0.20
round(cbind(noQTL_AC[,1],noQTL_nullcomp_AC),2)
#     [,1] [,2] [,3]
#[1,] 2.70 4.36 0.00
#[2,] 2.48 2.19 0.00
#[3,] 2.17 1.94 0.00
#[4,] 0.60 1.17 0.22
#[5,] 0.11 1.08 0.26
#[6,] 5.35 1.87 0.00
#[7,] 8.05 2.31 0.00

