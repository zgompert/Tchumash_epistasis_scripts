efGac<-scan("bmavAC.txt")

## genotypes accounting for centering
cGac<-matrix(c(-1.74489571,-1.87624877,-1.19373817,-1.27585320,-0.20343575,-0.05078484,
    -0.7448957,-0.8762488,-0.1937382,-0.2758532,0.7965642,0.9492152,
    0.2551043,0.1237512,0.8062618,0.7241468,1.7965642,1.9492152),nrow=3,byrow=TRUE)


blkEfAC<-c(0,-0.0309314931,-0.1884578527)
blkEfMM<-c(0,0.0064162945,0.0292963678)
b0<- mean(fit_AC[,3],na.rm=TRUE) + mean(blkEfAC)
#############################
# 2 way, try 2 x 5

gfit<-matrix(NA,nrow=3,ncol=3)
b2<-efGac[2]
b5<-efGac[5]
b2b5<-efGac[12]

for(i in 1:3){for(j in 1:3){
    gfit[i,j]<-b0 + b2 * cGac[i,2] + b5 * cGac[j,5] + b2b5 * cGac[i,2] * cGac[j,5]
}}

## ONLY BOTTOM 4 squares used
cs<-rev(heat.colors(n=9))
bs<-rev(c(1,0.7,.5,0.4,0.3,0.2,0.1,0,-0.3,-1))
pdf("epiExampleExfit2x5.pdf",width=4,height=4)
par(mar=c(4.5,4.5,0.5,0.5))
image(gfit,breaks=bs,col=cs,axes=FALSE,xlab="SNP 2",ylab="SNP 5")
        #title(main=paste("Fitness, PC1 = ",pcs[i],sep=""))
        axis(1,at=c(0,.5,1),c("aa","Aa","AA"))
        axis(2,at=c(0,.5,1),c("aa","Aa","AA"))
        box()

dev.off()
#plot(c(0,1),c(0,1),type='n',axes=F,xlab="",ylab="")
#image.plot(gfit[[i]],breaks=brks,col=cs,axes=FALSE,legend.only=TRUE,legend.lab="Rel. fitness")


# 3 way, try 1 x 3 x 4
gfit<-vector("list",3)
for(k in 1:3){
    gfit[[k]]<-matrix(NA,nrow=3,ncol=3)
}
b1<-efGac[1]
b3<-efGac[3]
b4<-efGac[4]
b1b3<-efGac[7]
b1b4<-efGac[8]
b3b4<-efGac[13]
b1b3b4<-efGac[19]

for(k in 1:3){
for(i in 1:3){for(j in 1:3){
    gfit[[k]][i,j]<-b0 + b1 * cGac[k,1] + b3 * cGac[i,3] + b4 * cGac[j,4] +
    b1b3 * cGac[k,1] * cGac[i,3] + b1b4 * cGac[k,1] * cGac[j,4] +
    b3b4 * cGac[i,3] * cGac[j,4] +
    b1b3b4 *  cGac[k,1] * cGac[i,3] * cGac[j,4]
  
}}
}

cs<-rev(heat.colors(n=9))
bs<-rev(c(1.5,0.7,.5,0.4,0.3,0.2,0.1,0,-0.3,-1))
pdf("epiExampleExfit1x3x4.pdf",width=8,height=8)
par(mfrow=c(2,2))
par(mar=c(4,4,2.5,2.5))
for(k in 1:3){
image(gfit[[k]],breaks=bs,col=cs,axes=FALSE,xlab="SNP 3",ylab="SNP 4")
        title(main=paste("SNP 1 = ",k-1,sep=""))
        axis(1,at=c(0,.5,1),c("aa","Aa","AA"))
        axis(2,at=c(0,.5,1),c("aa","Aa","AA"))
        box()
}
dev.off()
#plot(c(0,1),c(0,1),type='n',axes=F,xlab="",ylab="")
#image.plot(gfit[[i]],breaks=brks,col=cs,axes=FALSE,legend.only=TRUE,legend.lab="Rel. fitness")



# 3 way, color/fitness plot 
colEffs<-matrix(NA,nrow=3,ncol=2)
colEffs[1,]<- -1*c(1.023523e-05,2.136363e-02) # SNP1
colEffs[2,]<-c(-6.836257e-02,2.591891e-02) # SNP3
colEffs[3,]<-c(-3.198406e-03 -1.893820e-02) # SNP4
rg_cfit<-vector("list",3)
gb_cfit<-vector("list",3)
for(k in 1:3){
    rg_cfit[[k]]<-matrix(NA,nrow=3,ncol=3)
    gb_cfit[[k]]<-matrix(NA,nrow=3,ncol=3)
}

## for intercept on color
load("models.rdat")
cint<-matrix(NA,nrow=3,ncol=2)
cint[1,]<--1*c(-2.978231e-02,4.242169e-04) #SNP 2
cint[2,]<-c(4.027592e-02,-7.155026e-05) ## SNP 5
cint[3,]<-c(-6.979418e-02,1.055617e-01) ## SNP 6

p_ref<-(apply(g_AC_pip,2,mean)/2)[c(2,5,6)]
p_ref[1]<-1-p_ref[1]

## read in data
cdat<-read.table("../../color_exp_2019/2019_Tchumash_transplant_table.csv",header=TRUE,sep=",")

ac<-which(cdat$Treatment=='AC')

rg0<-mean(cdat$RG[ac]) + (2 * p_ref[1] * cint[1,1]) + (2 * p_ref[2] * cint[2,1]) +
	(2 * p_ref[3] * cint[3,1])
gb0<-mean(cdat$GB[ac]) + (2 * p_ref[1] * cint[1,2]) + (2 * p_ref[2] * cint[2,2]) +
	(2 * p_ref[3] * cint[3,2])

gfit<-vector("list",3)
for(k in 1:3){
    gfit[[k]]<-matrix(NA,nrow=3,ncol=3)
}
b1<-efGac[1]
b3<-efGac[3]
b4<-efGac[4]
b1b3<-efGac[7]
b1b4<-efGac[8]
b3b4<-efGac[13]
b1b3b4<-efGac[19]

genotypes<-vector("list",3)

gs<-c("aa","Aa","AA")
for(k in 1:3){
	genotypes[[k]]<-data.frame(matrix(NA,nrow=3,ncol=3))
for(i in 1:3){for(j in 1:3){
    gfit[[k]][i,j]<-b0 + b1 * cGac[k,1] + b3 * cGac[i,3] + b4 * cGac[j,4] +
    b1b3 * cGac[k,1] * cGac[i,3] + b1b4 * cGac[k,1] * cGac[j,4] +
    b3b4 * cGac[i,3] * cGac[j,4] +
    b1b3b4 *  cGac[k,1] * cGac[i,3] * cGac[j,4]

    genotypes[[k]][i,j]<-paste(c(gs[k],gs[i],gs[j]),collapse="")
    rg_cfit[[k]][i,j]<- (k-1) * colEffs[1,1] + (i-1) * colEffs[2,1] + (j-1) * colEffs[3,1]
    gb_cfit[[k]][i,j]<- (k-1) * colEffs[1,2] + (i-1) * colEffs[2,2] + (j-1) * colEffs[3,2]
    #rg_cfit[[k]][i,j]<- rg0 + (k-1) * colEffs[1,1] + (i-1) * colEffs[2,1] + (j-1) * colEffs[3,1]
    #gb_cfit[[k]][i,j]<- gb0+ (k-1) * colEffs[1,2] + (i-1) * colEffs[2,2] + (j-1) * colEffs[3,2]
  
}}
}

rg_vec<-c(as.vector(rg_cfit[[2]]),as.vector(rg_cfit[[3]]))
gb_vec<-c(as.vector(gb_cfit[[2]]),as.vector(gb_cfit[[3]]))
fit_vec<-b0+c(as.vector(gfit[[2]]),as.vector(gfit[[3]]))
g_vec<-c(as.vector(as.matrix(genotypes[[2]])),as.vector(as.matrix(genotypes[[3]])))

cs<-heat.colors(n=10)
cs_fit<-rep(cs[1],18)
brk<-c(1,0.65,.45,0.35,0.25,0.15,0.05,-0.05,-0.3)
for(i in 1:9){
    a<-which(fit_vec<brk[i])
    cs_fit[a]<-cs[i+1]
    }


pdf("epiScatter.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
plot(rg_vec,gb_vec,pch=20,col=cs_fit,xlab="Change in RG",ylab="Change in GB",cex.lab=1.5,
	xlim=c(-.2,0.02),cex=1.5)
text(rg_vec,gb_vec+0.003,g_vec,cex=.76)
dev.off()

cs<-rev(heat.colors(n=9))
bs<-rev(c(1.5,0.7,.5,0.4,0.3,0.2,0.1,0,-0.3,-1))
pdf("epiExampleExfit1x3x4.pdf",width=8,height=8)
par(mfrow=c(2,2))
par(mar=c(4,4,2.5,2.5))
for(k in 1:3){
image(gfit[[k]],breaks=bs,col=cs,axes=FALSE,xlab="SNP 3",ylab="SNP 4")
        title(main=paste("SNP 1 = ",k-1,sep=""))
        axis(1,at=c(0,.5,1),c("aa","Aa","AA"))
        axis(2,at=c(0,.5,1),c("aa","Aa","AA"))
        box()
}
dev.off()


