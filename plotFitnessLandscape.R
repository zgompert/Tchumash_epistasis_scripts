## plots the fitness landscape as a graph

library(igraph) 

# 3 binary
G<-cbind(rep(c(0,1),each=4),rep(c(0,0,1,1),2),rep(c(0,1),4))
nds<-paste("s",1:8,sep="")
colnames(nds)<-"id"
jns<-matrix(NA,nrow=28,ncol=3)
jns<-as.data.frame(jns)
n<-1
for(i in 1:7){for(j in (i+1):8){
    jns[n,1]<-nds[i]
    jns[n,2]<-nds[j]
    jns[n,3]<-sum(abs(G[i,]-G[j,]))
    n<-n+1
}}

colnames(jns)<-c("from","to","weight")
jns2<-jns[which(jns[,3]==1),]
net<-graph_from_data_frame(d=jns2,vertices=nds,directed=FALSE)
l<-layout_with_gem(net)
plot(net,layout=l)

# 5 binary
2^5
#[1] 32
G<-cbind(rep(c(0,1),each=16),
    rep(c(rep(0,8),rep(1,8)),2),
    rep(c(rep(0,4),rep(1,4)),4),
    rep(c(0,0,1,1),8),
    rep(c(0,1),16))

nds<-paste("s",1:32,sep="")
nds<-as.data.frame(nds)
colnames(nds)<-"id"
(32*31)/2
jns<-matrix(NA,nrow=496,ncol=3)
jns<-as.data.frame(jns)
n<-1
for(i in 1:31){for(j in (i+1):32){
    jns[n,1]<-as.character(nds[i,1])
    jns[n,2]<-as.character(nds[j,1])
    jns[n,3]<-sum(abs(G[i,]-G[j,]))
    n<-n+1
}}

colnames(jns)<-c("from","to","weight")
jns2<-jns[which(jns[,3]==1),]
net<-graph_from_data_frame(d=jns2,vertices=nds,directed=FALSE)
l<-layout_with_gem(net)
plot(net,layout=l)

# 5, 3 genotype
3^5
#[1] 243
G<-cbind(rep(c(0,1,2),each=81),
    rep(c(rep(0,27),rep(1,27),rep(2,27)),3),
    rep(c(rep(0,9),rep(1,9),rep(2,9)),9),
    rep(c(0,0,0,1,1,1,2,2,2),27),
    rep(c(0,1,2),81))

nds<-rep(NA,243)
for(i in 1:243){
    nds[i]<-paste(G[i,],collapse="")
    }
nds<-as.data.frame(nds)
colnames(nds)<-"id"
(243*242)/2
jns<-matrix(NA,nrow=29403,ncol=3)
jns<-as.data.frame(jns)
n<-1
for(i in 1:242){for(j in (i+1):243){
    jns[n,1]<-as.character(nds[i,1])
    jns[n,2]<-as.character(nds[j,1])
    jns[n,3]<-sum(abs(G[i,]-G[j,]))
    n<-n+1
}}

colnames(jns)<-c("from","to","weight")
jns2<-jns[which(jns[,3]==1),]
net<-graph_from_data_frame(d=jns2,vertices=nds,directed=FALSE)
l<-layout_with_gem(net)
plot(net,layout=l)
plot(net,layout=l,vertex.size=2,vertex.color=c(rep("orange",100),rep("yellow",143)))

## model averaged effects from all epistasis model
efGac<-scan("bmavAC.txt")

## genotypes accounting for centering
cGac<-matrix(c(-1.74489571,-1.87624877,-1.19373817,-1.27585320,-0.20343575,-0.05078484,
    -0.7448957,-0.8762488,-0.1937382,-0.2758532,0.7965642,0.9492152,
    0.2551043,0.1237512,0.8062618,0.7241468,1.7965642,1.9492152),nrow=3,byrow=TRUE)

exfit<-rep(NA,243)
for(i in 1:243){
    gg<-rep(NA,5)
    for(j in 1:5){
        gg[j]<-cGac[G[i,j]+1,j]
    }
    y<-1
    form<-y~.^5
    mod<-model.matrix(form,dat=as.data.frame(t(gg)))
    exfit[i]<-sum(mod[-1] * efGac)
    }

exfitAC<-exfit

## model averaged effects from all epistasis model
efGmm<-scan("bmavMM.txt")

## genotypes accounting for centering
cGmm<-matrix(c(-1.74272546,-1.83116023,-1.20107060,-1.38443852,-0.15445727,-0.04479259,
    -0.7427255,-0.8311602,-0.2010706,-0.3844385,0.8455427,0.9552074,
    0.2572745,0.1688398,0.7989294,0.6155615,1.8455427,1.9552074),nrow=3,byrow=TRUE)

exfit<-rep(NA,243)
for(i in 1:243){
    gg<-rep(NA,5)
    for(j in 1:5){
        gg[j]<-cGmm[G[i,j]+1,j]
    }
    y<-1
    form<-y~.^5
    mod<-model.matrix(form,dat=as.data.frame(t(gg)))
    exfit[i]<-sum(mod[-1] * efGmm)
    }

exfitMM<-exfit

## count obs.
rg_AC_pip<-round(g_AC_pip)
N<-dim(rg_AC_pip)[1]
cntsAC<-rep(0,243)
for(i in 1:N){
	xm<-which(nds==paste(rg_AC_pip[i,-6],collapse=""))
	cntsAC[xm]<-1 + cntsAC[xm]
}
pchAC<-rep(2,243)
pchAC[cntsAC > 0]<-4.5
pchAC[cntsAC > 5]<-6

rg_MM_pip<-round(g_MM_pip)
N<-dim(rg_MM_pip)[1]
cntsMM<-rep(0,243)
for(i in 1:N){
	xm<-which(nds==paste(rg_MM_pip[i,-6],collapse=""))
	cntsMM[xm]<-1 + cntsMM[xm]
}
pchMM<-rep(2,243)
pchMM[cntsMM > 0]<-4.5
pchMM[cntsMM > 5]<-6

library(RColorBrewer)    

load("models.rdat")

blkEfAC<-c(0,-0.0309314931,-0.1884578527)
blkEfMM<-c(0,0.0064162945,0.0292963678)

exfitAC<-exfitAC + mean(fit_AC[,3],na.rm=TRUE) + mean(blkEfAC)
exfitMM<-exfitMM + mean(fit_MM[,3],na.rm=TRUE) + mean(blkEfMM)

pdf("chumashFitLands.pdf",width=22,height=11)
cs<-heat.colors(n=10)
#cs<-rev(brewer.pal(n=9,"RdYlGn"))
par(mfrow=c(1,2))
csfit<-rep(cs[1],243)
brk<-c(1,0.65,.45,0.35,0.25,0.15,0.05,-0.05,-0.3)
#brk<-c(1,0.7,.5,0.4,0.3,0.2,0.1,0,-0.3,-0.5)
for(i in 1:9){
    a<-which(exfitAC<brk[i])
    csfit[a]<-cs[i+1]
    }
plot(net,layout=l,vertex.label=NA,vertex.size=pchAC,vertex.color=csfit)
title(main="(A) A/C fitness landscape",cex.main=1.7)

csfit<-rep(cs[1],243)
brk<-c(1,0.65,.45,0.35,0.25,0.15,0.05,-0.05,-0.3)
#brk<-c(1,0.7,.5,0.4,0.3,0.2,0.1,0,-0.3,-0.5)
for(i in 1:9){
    a<-which(exfitMM<brk[i])
    csfit[a]<-cs[i+1]
    }
plot(net,layout=l,vertex.label=NA,vertex.size=pchMM,vertex.color=csfit)
title(main="(B) MM fitness landscape",cex.main=1.7)
dev.off()

## NOTE fitness values need means added if you want them to be absolute,
## and then they are really for block 1, need additional terms for blocks 2 and 3

## count peaks
pkComp<-rep(NA,243)
for(i in 1:243){
	Gref<-G[i,]
		Yf<-exfitAC[i]
	Yxf<-NULL
	for(k in 1:5){
		if((Gref[k] == 0) | (Gref[k]==2)){
			GG<-Gref
			GG[k]<-1
			GGp<-paste(GG,collapse="")
			Yxf<-c(Yxf,exfitAC[as.character(nds[,1])==as.character(GGp)])
		}
		else if(Gref[k] == 1){
			GG<-Gref
			GG[k]<-0
			GGp<-paste(GG,collapse="")
			Yxf<-c(Yxf,exfitAC[as.character(nds[,1])==as.character(GGp)])
			GG[k]<-2
			GGp<-paste(GG,collapse="")
			Yxf<-c(Yxf,exfitAC[as.character(nds[,1])==as.character(GGp)])
		}
	}
	pkComp[i]<-mean(Yxf < Yf)
}
sum(pkComp==1)

pkComp<-rep(NA,243)
for(i in 1:243){
	Gref<-G[i,]
	Yf<-exfitMM[i]
	Yxf<-NULL
	for(k in 1:5){
		if((Gref[k] == 0) | (Gref[k]==2)){
			GG<-Gref
			GG[k]<-1
			GGp<-paste(GG,collapse="")
			Yxf<-c(Yxf,exfitMM[as.character(nds[,1])==as.character(GGp)])
		}
		else if(Gref[k] == 1){
			GG<-Gref
			GG[k]<-0
			GGp<-paste(GG,collapse="")
			Yxf<-c(Yxf,exfitMM[as.character(nds[,1])==as.character(GGp)])
			GG[k]<-2
			GGp<-paste(GG,collapse="")
			Yxf<-c(Yxf,exfitMM[as.character(nds[,1])==as.character(GGp)])
		}
	}
	pkComp[i]<-mean(Yxf < Yf)
}
sum(pkComp==1)

## variances
sd(exfitAC[cntsAC>0])
#[1] 0.187152
sd(exfitMM[cntsMM>0])
#[1] 0.04051042
sd(exfitAC)
#[1] 0.691771
sd(exfitMM)
#[1] 0.2857427

## random walks for ruggedness
rwalk<-function(GG=G,start=1,nstep=10,nds=nds){
	pos<-rep(NA,nstep+1)
	pos[1]<-start
	for(i in 1:nstep){
		mv<-FALSE
		while(mv==FALSE){
			k<-sample(1:5,1)
			d<-sample(c(-1,1),1)
			Gx<-G[pos[i],]
			Gx[k]<-G[pos[i],k]+d
			if(sum(Gx >= 0 & Gx <= 2)==5){
				mv<-TRUE
			}
		}
		pos[i+1]<-which(as.character(nds[,1])==as.character(paste(Gx,collapse="")))
	}
	return(pos)
}


oc<-which(cntsAC > 0 | cntsMM > 0)

## 10 walk, from occ.
scAC<-rep(NA,10000)
scMM<-rep(NA,10000)
rscAC<-rep(NA,10000)
rscMM<-rep(NA,10000)
for(x in 1:10000){
	st<-sample(oc,1)
	o<-exfitAC[rwalk(GG=G,start=st,nstep=10,nds=nds)]
	dif<-o[-1]-o[-10]
	scAC[x]<-sum(abs(dif))-abs(sum(dif))
	rscAC[x]<-sum(abs(dif))/abs(sum(dif))
	o<-exfitMM[rwalk(GG=G,start=st,nstep=10,nds=nds)]
	dif<-o[-1]-o[-10]
	scMM[x]<-sum(abs(dif))-abs(sum(dif))
	rscMM[x]<-sum(abs(dif))/abs(sum(dif))
	}

summary(scAC);sd(scAC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.06884  1.43183  2.04446  2.25731  2.81465 10.71626 
#[1] 1.203165
summary(rscAC);sd(rscAC)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#    1.05     4.02     7.49    64.08    16.42 52008.06 
#[1] 1015.606
summary(scMM);sd(scMM)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.005935 0.159728 0.314929 0.494732 0.619255 6.031544 
#[1] 0.546643
summary(rscMM);sd(rscMM)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   1.038    3.911    7.969   33.276   19.246 8372.521 
#[1] 170.5005


## 10 walk, from randome
scAC<-rep(NA,10000)
scMM<-rep(NA,10000)
rscAC<-rep(NA,10000)
rscMM<-rep(NA,10000)
for(x in 1:10000){
	st<-sample(1:243,1)
	o<-exfitAC[rwalk(GG=G,start=st,nstep=10,nds=nds)]
	dif<-o[-1]-o[-10]
	scAC[x]<-sum(abs(dif))-abs(sum(dif))
	rscAC[x]<-sum(abs(dif))/abs(sum(dif))
	o<-exfitMM[rwalk(GG=G,start=st,nstep=10,nds=nds)]
	dif<-o[-1]-o[-10]
	scMM[x]<-sum(abs(dif))-abs(sum(dif))
	rscMM[x]<-sum(abs(dif))/abs(sum(dif))
}
summary(scAC);sd(scAC)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.03707  1.84275  2.75332  3.06286  3.93859 13.27581 
#[1] 1.686596
summary(rscAC);sd(rscAC)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#    1.03     3.49     6.24    39.11    13.82 46010.26 
#[1] 577.7476
summary(scMM);sd(scMM)
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.001745  0.331231  0.702296  1.119944  1.465065 11.291032 
#[1] 1.19294
summary(rscMM);sd(rscMM)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#    1.01     3.81     7.61    42.29    17.33 42852.47 
#[1] 627.2249

summary(scAC);sd(scAC)
summary(rscAC);sd(rscAC)
summary(scMM);sd(scMM)
summary(rscMM);sd(rscMM)

## 20 walk, from occ.
scAC<-rep(NA,10000)
scMM<-rep(NA,10000)
rscAC<-rep(NA,10000)
rscMM<-rep(NA,10000)
for(x in 1:10000){
	st<-sample(oc,1)
	o<-exfitAC[rwalk(GG=G,start=st,nstep=20,nds=nds)]
	dif<-o[-1]-o[-20]
	scAC[x]<-sum(abs(dif))-abs(sum(dif))
	rscAC[x]<-sum(abs(dif))/abs(sum(dif))
	o<-exfitMM[rwalk(GG=G,start=st,nstep=20,nds=nds)]
	dif<-o[-1]-o[-20]
	scMM[x]<-sum(abs(dif))-abs(sum(dif))
	rscMM[x]<-sum(abs(dif))/abs(sum(dif))
	}

summary(scAC);sd(scAC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7817  4.4225  5.7454  6.1076  7.3913 20.7864 
#[1] 2.346485
summary(rscAC);sd(rscAC)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#     1.74      9.11     17.07    183.96     38.52 168210.43 
#[1] 3349.834
summary(scMM);sd(scMM)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0604  0.7114  1.2608  1.7223  2.1988 13.4518 
#[1] 1.518863
summary(rscMM);sd(rscMM)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#    1.265     9.642    21.881   114.420    55.160 19661.150 
#[1] 579.4403

## Moran I
N<-243
#AC
ybar<-mean(exfitAC)
dy<-(exfitAC-ybar)
g<-expand.grid(dy,dy)
yiyj<-g[,1]*g[,2]
pm<-matrix(yiyj,N)
w<-matrix(0,nrow=N,ncol=N)
for(i in 1:243){for(j in 1:243){
	if(abs(sum(G[i,]-G[j,]))==1){
		w[i,j]<-1
	}
}}
pmw<-pm*w
spmw<-sum(pmw)
sw<-spmw/(sum(w))
vr<-N/sum(dy^2)
MI_AC<-vr*sw
## subset with genotypes
xpmw<-pmw[cntsAC>0,cntsAC>0]
sxpmw<-sum(xpmw)
xsw<-sxpmw/sum(w[cntsAC>0,cntsAC>0])
xvr<-sum(cntsAC>0)/sum(dy[cntsAC>0]^2)
xMI_AC<-xvr*xsw
xMI_AC
#[1] 0.09799135
MI_AC
#[1] 0.02875857


##MM
ybar<-mean(exfitMM)
dy<-(exfitMM-ybar)
g<-expand.grid(dy,dy)
yiyj<-g[,1]*g[,2]
pm<-matrix(yiyj,N)
w<-matrix(0,nrow=N,ncol=N)
for(i in 1:243){for(j in 1:243){
	if(abs(sum(G[i,]-G[j,]))==1){
		w[i,j]<-1
	}
}}
pmw<-pm*w
spmw<-sum(pmw)
sw<-spmw/(sum(w))
vr<-N/sum(dy^2)
MI_MM<-vr*sw
## subset with genotypes
xpmw<-pmw[cntsMM>0,cntsMM>0]
sxpmw<-sum(xpmw)
xsw<-sxpmw/sum(w[cntsMM>0,cntsMM>0])
xvr<-sum(cntsMM>0)/sum(dy[cntsMM>0]^2)
xMI_MM<-xvr*xsw

MI_MM
#[1] 0.03201357
xMI_MM
#[1] 0.06100262
