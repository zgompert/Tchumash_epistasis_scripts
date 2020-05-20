library(scales)
library(latex2exp)

rgPip<-scan("output/pip_o_mod_g_tchum_ph1_ch0.param.txt")
gbPip<-scan("output/pip_o_mod_g_tchum_ph2_ch0.param.txt")
acPip<-scan("output/pip_o_mod_g_tchum_AC_ph2_ch0.param.txt")
mmPip<-scan("output/pip_o_mod_g_tchum_MM_ph2_ch0.param.txt")


snps<-read.table("snpIds.txt",header=FALSE)
ss<-read.table("snpSorter.txt",header=FALSE)
snps<-cbind(snps,rep(NA,dim(snps)[1]))
for(i in 1:dim(snps)[1]){
	a<-which(ss[,3]==snps[i,2])
	if(length(a)==1){
		snps[i,4]<-ss[a,2]
	}
}


###################
# compute pp for no. of QTN
# focus on mel-stripe
ms<-which((snps[,2] == 702.1 & snps[,3] <= 4139489) | (snps[,2] == 128 & snps[,3] < 6414835))

qtnRG<-rep(NA,10000)
qtnGB<-rep(NA,10000)
qtnCom<-rep(NA,10000)
Nms<-length(ms)
for(i in 1:10000){
	rg<-rbinom(n=Nms,size=1,prob=rgPip[ms])	
	gb<-rbinom(n=Nms,size=1,prob=gbPip[ms])
	rbgb<-apply(cbind(rg,gb),1,max)
	qtnRG[i]<-sum(rg)	
	qtnGB[i]<-sum(gb)
	qtnCom[i]<-sum(rbgb)
}	
mean(qtnRG)
#[1] 4.6293
sd(qtnRG)
#[1] 0.8082987
mean(qtnGB)
#[1] 3.9406
sd(qtnGB)
#[1] 1.034784
mean(qtnCom)
#[1] 6.5079
sd(qtnCom)
#[1] 1.152333


###################
mkplotInDel<-function(X=NA,c1=NA,snp=snps,hdr=NA){
	indel<-which(snp[,2] == 128 & snp[,3] >= 5000000 & snp[,3] <= 6000000)

	cx<-1.6
	ca<-1.2
	cm<-1.6
	mypc<-rep(20,length(X))
	mypc[X>.5]<-19
	plot(snp[indel,3],X[indel],col=c1,pch=mypc[indel],axes=FALSE,ylim=c(0,1),cex=1.1,cex.lab=cx,xlab="Scaffold 128 pos. (bp)",ylab="PIP")
	axis(1,cex.axis=ca)
	axis(2,at=c(0,.5,1),cex.axis=ca)
	title(main=hdr,cex.main=cm)
	box()
}
estPipInDel<-function(X=NA,snp=snps){
	
	indel<-which(snp[,2] == 128 & snp[,3] >= 5000000 & snp[,3] <= 6000000)
	Ni<-length(indel)
	Nq<-rep(NA,10000)
	for(k in 1:10000){
		Nq[k]<-sum(rbinom(n=Ni,size=1,prob=X[indel]))
	}
	o<-c(mean(Nq),sd(Nq))
	return(o)
}


mkplotGw<-function(X=NA,c2=NA,hdr=NA,snp=snps){
	drop<-which(snp[,1]==0)
	X<-X[-drop]
	snp<-snp[-drop,]
	X<-X[order(snp[,1],snp[,4])]
	snp<-snp[order(snp[,1],snp[,4]),]
	L<-length(X)
	minp<-which(X>0.0005)
	ps<-1:L
	mids<-tapply(X=ps,INDEX=snp[,1],mean)
	cx<-1.6
	ca<-1.2
	cm<-1.6
	mycs<-rep(c2[1],L)
	mycs[snp[,1]==8]<-c2[2]
	mypc<-rep(20,L)
	mypc[X>.5]<-19
	plot(ps[minp],X[minp],col=mycs[minp],pch=mypc[minp],axes=FALSE,ylim=c(0,1),cex.lab=cx,xlab="Linkage group",ylab="PIP")
	axis(1,at=mids,1:13,cex.axis=ca)
	axis(2,at=c(0,.5,1),cex.axis=ca)
	title(main=hdr,cex.main=cm)
	box()
}
	

## make plots
pdf("chumashColorPips.pdf",width=11,height=10)
par(mfrow=c(2,2))

cs<-c("darkgray","indianred2","cornflowerblue")

par(mar=c(4.5,5.5,3,.5))

## pip plots
mkplotGw(X=rgPip,c2=cs[1:2],hdr=expression(paste("(a) Genome-wide PIPs, RG")))
mkplotInDel(X=rgPip,c1=cs[2],hdr=expression(paste("(b) Scaf. 128 PIPs, RG")))
mkplotGw(X=gbPip,c2=cs[c(1,3)],hdr=expression(paste("(c) Genome-wide PIPs, GB")))
mkplotInDel(X=gbPip,c1=cs[c(3)],hdr=expression(paste("(d) Scaf. 128 PIPs, GB")))
dev.off()

## make plots
pdf("chumashSurvPips.pdf",width=11,height=10)
par(mfrow=c(2,2))

cs<-c("darkgray","forestgreen","forestgreen")

par(mar=c(4.5,5.5,3,.5))

## pip plots
mkplotGw(X=acPip,c2=cs[1:2],hdr=expression(paste("(a) Genome-wide PIPs, A/C")))
mkplotInDel(X=acPip,c1=cs[2],hdr=expression(paste("(b) Scaf. 128 PIPs, A/C")))
mkplotGw(X=mmPip,c2=cs[c(1,3)],hdr=expression(paste("(c) Genome-wide PIPs, MM")))
mkplotInDel(X=mmPip,c1=cs[c(3)],hdr=expression(paste("(d) Scaf. 128 PIPs, MM")))
dev.off()
