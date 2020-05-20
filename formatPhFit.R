## line up expected fitness with phenotype/genotype data

N<-454

keep<-scan("../Variants/hiCov2x.txt") ## single column, denotes whether 1 or not 0 mean cov for an individual > 2
## true for 437, going to drop the others

gids<-read.table("../Variants/IdsGen.txt",header=FALSE)

## read phenotype and fitness data
phtab<-read.table("../Variants/2019_Tchumash_transplant_table.csv",sep=",",header=TRUE)
phtab2<-read.table("../../color_exp_2019/2019_Tchumash_transplant_table.csv",sep=",",header=TRUE)
fit<-read.table("fitAll.txt",header=FALSE)

ph<-matrix(NA,nrow=N,ncol=8)
ph<-as.data.frame(ph)
for(i in 1:N){
	a<-which(as.character(phtab[,1])==as.character(gids[i,1]))
	if(length(a)==1){
		ph[i,]<-phtab[a,c(1,3:6,11:13)]
	}
}
colnames(ph)<-colnames(phtab)[c(1,3:6,11:13)]
#number codes
# Geneder: M = 2, F = 1
# Color: G = 1, GB = 2, I = 3, M = 4
# Treatment: AC = 1, MM = 2

## fit
fitt<-matrix(NA,nrow=N,ncol=4)
fitt<-as.data.frame(fitt)
for(i in 1:N){
	a<-which(as.character(phtab2[,1])==as.character(gids[i,1]))
	if(length(a)==1){
		fitt[i,]<-fit[a,]
	}
}

## subset those that have good genetic data
g_hc<-g[,keep==1]
ph_hc<-ph[keep==1,]
fitt_hc<-fitt[keep==1,]
dim(g_hc)
#[1] 11990   437
dim(ph_hc)
#[1] 437   8
dim(fitt_hc)
#[1] 437   4
fitt_hc[,4]<-ph_hc$Treatment

## now by treatment, gender and survival (former is kind of for fun)
AC<-which(fitt_hc[,4]==1)
MM<-which(fitt_hc[,4]==2)

colnames(fitt_hc)<-c("survival","blk","exfit","tretament")
write.table(file="exfit_AC.txt",fitt_hc[AC,1:3],row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(file="exfit_MM.txt",fitt_hc[MM,1:3],row.names=FALSE,col.names=TRUE,quote=FALSE)

