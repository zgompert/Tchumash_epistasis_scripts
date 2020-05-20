## create effect size plots from sampling posterior 

b_rg<-scan("output/rawbeta_o_mod_g_tchum_ph1_ch0.param.txt")
b_gb<-scan("output/rawbeta_o_mod_g_tchum_ph2_ch0.param.txt")
pip_rg<-scan("output/pip_o_mod_g_tchum_ph1_ch0.param.txt")
pip_gb<-scan("output/pip_o_mod_g_tchum_ph2_ch0.param.txt")

L<-11855
N<-10000

efd_rg<-NULL
efd_gb<-NULL

## genome wide

for(i in 1:N){
	efd_rg<-c(efd_rg,b_rg[rbinom(n=L,size=1,prob=pip_rg)==1])
	efd_gb<-c(efd_gb,b_gb[rbinom(n=L,size=1,prob=pip_gb)==1])
}

pdf("chumColorEffDist.pdf",width=9,height=5)
cl<-1.5
cm<-1.3
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,3,1))
cs<-c("indianred2","cornflowerblue")
hist(abs(efd_rg),col=cs[1],xlab="Effect size",ylab="Frequency",main="(A) Eff. size distr. RG",cex.lab=cl,cex.main=cm,xlim=c(0,0.12))
box()
hist(abs(efd_gb),col=cs[2],xlab="Effect size",ylab="Frequency",main="(B) Eff. size distr. GB",cex.lab=cl,cex.main=cm,xlim=c(0,0.12))
box()
dev.off()

## mel stripe only
snps<-read.table("output/pipSnpTable.txt",header=FALSE)

ms<-c(which(snps[,2]==128 & snps[,3] <6000000),which(snps[,2]==702.1 & snps[,3] < 4139489))
Lms<-157


N<-10000

efd_rg<-NULL
efd_gb<-NULL

## genome wide

for(i in 1:N){
	efd_rg<-c(efd_rg,b_rg[ms][rbinom(n=Lms,size=1,prob=pip_rg[ms])==1])
	efd_gb<-c(efd_gb,b_gb[ms][rbinom(n=Lms,size=1,prob=pip_gb[ms])==1])
}

pdf("chumColorEffDistMs.pdf",width=9,height=5)
cl<-1.5
cm<-1.3
par(mfrow=c(1,2))
par(mar=c(4.5,4.5,3,1))
cs<-c("indianred2","cornflowerblue")
hist(abs(efd_rg),col=cs[1],xlab="Effect size",ylab="Frequency",main="(A) Eff. size distr. RG",cex.lab=cl,cex.main=cm,xlim=c(0,0.12))
box()
hist(abs(efd_gb),col=cs[2],xlab="Effect size",ylab="Frequency",main="(B) Eff. size distr. GB",cex.lab=cl,cex.main=cm,xlim=c(0,0.12))
box()
dev.off()
