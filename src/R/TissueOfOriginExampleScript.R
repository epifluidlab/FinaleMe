# TissueOfOriginExampleScript.R
# Aug 17, 2016
# 4:32:54 PM
# 
# Author: yaping
###############################################################################

library(quadprog)
tissue_of_origin<-function(dat){
	dat<-dat[rowSums(is.na(dat))==0 & rowSums(is.infinite(as.matrix(dat)))==0,]

	X<-as.matrix(cbind(dat[,2:length(dat[1,])]))
	Y<-dat[,1]

	Rinv <- solve(chol(t(X) %*% X));
	C <- cbind(rep(1,dim(X)[2]), diag(dim(X)[2]))
	b <- c(1,rep(0,dim(X)[2]))
	d <- t(Y) %*% X
	result<-solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec=b, meq=1)
	result$solution
}

setwd("/jet/home/dnaase/startup/projects/ccinference/broad_data/wgs_bam/too")
raw.1kb<-read.table("b37.autosome.cpgIsland_plus_shore_1kb_intervals.add_value.methy.wgbs_ref_14_samples.bed.gz",sep="\t",header=F)
name.order<-read.table("names_order.ref_panel_breast_prostate.txt",sep="\t",header=F)

methy.mat<-NULL
for(i in seq(7,length(raw.1kb[1,]),2)){
	j=i+1
	cov = raw.1kb[,j]
	cov[cov*1000<10]=0
	s=as.numeric(raw.1kb[,i]/cov)
	methy.mat<-cbind(methy.mat,s)
}
rownames(methy.mat)<-gsub("chr","",raw.1kb[,4])
colnames(methy.mat)<-name.order[,2]

methy.mat=methy.mat[,-c(2,4,5,13)]
name.order=name.order[-c(2,4,5,13),]

row.sd<-apply(methy.mat,1, sd, na.rm = TRUE)
quantile(row.sd,seq(0,1,0.1), na.rm=T)

methy.mat.mostVar=methy.mat[row.sd>quantile(row.sd,seq(0,1,0.01), na.rm=T)[100],]
methy.mat.mostVar=methy.mat.mostVar[rowSums(is.na(methy.mat.mostVar),na.rm=T)==0,]
methy.mat.mostVar=methy.mat.mostVar[rowSums(is.infinite(methy.mat.mostVar),na.rm=T)==0,]
methy.mat.mostVar[methy.mat.mostVar<0.1]=0
methy.mat.mostVar[methy.mat.mostVar>=0.1]=1

methy.mat.mostVar=methy.mat.mostVar[,c(6,1,7,5,4,2,8,3,9,10)]
name.order=name.order[c(6,1,7,5,4,2,8,3,9,10),]

cfdna.1kb<-read.table("b37.autosome.cpgIsland_plus_shore_1kb_intervals.self_train_mincg10.bin1.add_value.methy.broad_deep_wgs_wgbs_trim_3end_50bp.bed.gz",sep="\t",header=F)
cfdna.name.order<-read.table("names_order.deep_wgs_wgbs.txt",sep="\t",header=F)

methy.cfmat<-NULL
for(i in seq(7,length(cfdna.1kb[1,]),2)){
	j=i+1
	cov = cfdna.1kb[,j]
	cov[cov*1000<10]=0
	s=as.numeric(cfdna.1kb[,i]/cov)
	methy.cfmat<-cbind(methy.cfmat,s)
}
rownames(methy.cfmat)<-cfdna.1kb[,4]
colnames(methy.cfmat)<-cfdna.name.order[,2]
methy.cfmat=methy.cfmat[rowSums(is.na(methy.cfmat),na.rm=T)==0,]
methy.cfmat=methy.cfmat[rowSums(is.infinite(methy.cfmat),na.rm=T)==0,]
methy.cfmat[methy.cfmat<0.1]=0
methy.cfmat[methy.cfmat>=0.1]=1

common<-intersect(rownames(methy.cfmat),rownames(methy.mat.mostVar))
methy.cfmat.common<-methy.cfmat[common,]
methy.mat.ref<-methy.mat.mostVar[common,]

too_result.1kb<-NULL
for(i in 1:length(methy.cfmat.common[1,])){
	dat<-cbind(methy.cfmat.common[,i],methy.mat.ref)
	res<-tissue_of_origin(dat)
	too_result.1kb<-rbind(too_result.1kb,res)
}
colnames(too_result.1kb)<-name.order[,2]
rownames(too_result.1kb)<-cfdna.name.order[,2]
too_result.1kb[too_result.1kb<0.001]=0