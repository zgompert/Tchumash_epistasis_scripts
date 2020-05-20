## set up phenotype files for 10 fold cross-validation

ph<-read.table("pheno_color.txt",header=FALSE)

rg_cv<-matrix(ph[,1],nrow=437,ncol=10)
gb_cv<-matrix(ph[,2],nrow=437,ncol=10)

grp<-sample(rep(1:10,each=44),437,replace=FALSE)

for(i in 1:10){
	xna<-which(grp==i)
	rg_cv[xna,i]<-NA
	gb_cv[xna,i]<-NA
}

write.table(rg_cv,file="pheno_rg_cv.txt",row.names=FALSE,col.names=FALSE)
write.table(gb_cv,file="pheno_gb_cv.txt",row.names=FALSE,col.names=FALSE)

