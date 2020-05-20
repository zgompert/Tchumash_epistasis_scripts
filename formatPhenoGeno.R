## format genotype data and phenotype data for gemma

## read in genotype point estimates
L<-11990
N<-454
g<-matrix(scan("pntest_filtered2x_tchum_mel_gbs_2019.txt",n=N*L,sep=" "),nrow=L,ncol=N,byrow=TRUE)

keep<-scan("../Variants/hiCov2x.txt") ## single column, denotes whether 1 or not 0 mean cov for an individual > 2
## true for 437, going to drop the others

gids<-read.table("../Variants/IdsGen.txt",header=FALSE)

## read phenotype data
phtab<-read.table("../Variants/2019_Tchumash_transplant_table.csv",sep=",",header=TRUE)
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

## subset those that have good genetic data
g_hc<-g[,keep==1]
ph_hc<-ph[keep==1,]

dim(g_hc)
#[1] 11990   437
dim(ph_hc)
#[1] 437   8

## write output for color mapping
## no need to control for gender, explains 1.8% of RG and <1% for GB
write.table(file="pheno_color.txt",ph_hc[,7:8],row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(file="g_tchum.txt",g_hc,row.names=FALSE,col.names=FALSE,quote=FALSE)

## now by treatment, gender and survival (former is kind of for fun)
AC<-which(ph_hc$Treatment==1)
MM<-which(ph_hc$Treatment==2)

write.table(file="pheno_AC.txt",ph_hc[AC,c(3,6)],row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(file="g_tchum_AC.txt",g_hc[,AC],row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table(file="pheno_MM.txt",ph_hc[MM,c(3,6)],row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(file="g_tchum_MM.txt",g_hc[,MM],row.names=FALSE,col.names=FALSE,quote=FALSE)

## note, survival is reasonably balanced
apply(ph_hc[MM,c(3,6)],2,mean)
apply(ph_hc[AC,c(3,6)],2,mean)
