# TissueOfOriginExampleScript.R
# Aug 17, 2016
# 4:32:54 PM
# 
# Author: yaping
###############################################################################


percBestGlm<-function(obs,ref){
	dataF<-data.frame(ref=ref,obs=obs)
	perc.predict=rep(0,dim(ref)[2])
	names(perc.predict)=colnames(dataF)[1:dim(ref)[2]]
	
	ans<-bestglm(dataF, IC="BICq")
	selectNames<-names(ans$BestModel$coefficients)
	selectNames<-selectNames[2:length(selectNames)]
	
	X <- as.matrix(dataF[,selectNames])
	Y <- as.matrix(obs)
	Rinv <- solve(chol(t(X) %*% X));
	C <- cbind(rep(1,length(selectNames)), diag(length(selectNames)))
	b <- c(1,rep(0,length(selectNames)))
	d <- t(Y) %*% X  
	result<-solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
	perc.predict[selectNames]=as.numeric(result$solution)
	perc.predict
}

setwd("/home/unix/li/compbio/project/cfDNA/simulation_mixture/encode_mix/mix_3_samples")
library(quadprog)
s1<-read.table("GM12878.ENCFF681ASN.GRCh38.mean_methy.no_overlap.window-100000.step-100000.min_data-100.bed",sep="\t",header=F)
rownames(s1)<-s1[,4]
s2<-read.table("K562.ENCFF963XLT.GRCh38.mean_methy.no_overlap.window-100000.step-100000.min_data-100.bed",sep="\t",header=F)
rownames(s2)<-s2[,4]
s3<-read.table("HepG2.ENCFF957OIM.GRCh38.mean_methy.no_overlap.window-100000.step-100000.min_data-100.bed",sep="\t",header=F)
rownames(s3)<-s3[,4]
s4<-read.table("H1.ENCFF546TLK.GRCh38.mean_methy.no_overlap.window-100000.step-100000.min_data-100.bed",sep="\t",header=F)
rownames(s4)<-s4[,4]

common<-intersect(rownames(s1),rownames(s2))
common<-intersect(common,rownames(s3))
common<-intersect(common,rownames(s4))

ref=data.frame(s1=s1[common,5], s2=s2[common,5],s3=s3[common,5],s4=s4[common,5])
rownames(ref)=common

predict_ratio.no_overlap<-NULL
obs_ratio.no_overlap<-NULL
rmses.no_overlap<-NULL
for(f in list.files(path = ".", pattern = "sampled_3000000.window-100000.step-100000.min_data-100.bed")){
	data<-read.table(f,sep="\t",header=F)
	rownames(data)<-data[,4]
	common.file<-intersect(common,rownames(data))
	obs<-data[common.file,5]
	result<-percBestGlm(ref=ref[common.file,],obs=obs)
	predict_ratio.no_overlap<-rbind(predict_ratio.no_overlap,as.numeric(result))
	
	r1<-as.numeric(gsub('GM12878.ENCFF681ASN.(.*).K562.ENCFF963XLT.*','\\1',f))
	r2<-as.numeric(gsub('.*.K562.ENCFF963XLT.(.*).HepG2.ENCFF957OIM.*','\\1',f))
	if(is.na(r2)){
		r2<-as.numeric(gsub('.*.K562.ENCFF963XLT.(.*).H1.ENCFF546TLK.*','\\1',f))
	}
	r3<-as.numeric(gsub('.*.HepG2.ENCFF957OIM.(.*).sampled_3000000.*','\\1',f))
	if(is.na(r3)){
		r3=0
	}
	r4<-as.numeric(gsub('.*.H1.ENCFF546TLK.(.*).sampled_3000000.*','\\1',f))
	if(is.na(r4)){
		r4=0
	}
	rmses.no_overlap<-c(rmses.no_overlap,sqrt(((r1-result[1])^2+(r2-result[2])^2 + (r3-result[3])^2 + (r4-result[4])^2)/4))
	obs_ratio.no_overlap<-rbind(obs_ratio.no_overlap,c(r1,r2,r3,r4))
	print(f)
	print(result)
}


pdf("tissue_of_origin_20160817.heatmap.4tissues_encode.infer_best_subset.predicted.pdf", paper="special", height=5, width=20)
par(oma=c(0, 0, 0, 0),mar=c(5, 5, 3, 1))
m<-t(apply(predict_ratio.no_overlap, 1, function(x)abs(x)/sum(abs(x),na.rm=T)))
image(1:nrow(m), 1:ncol(m), m, col = colorRampPalette(c("white","red"))(100), axes = FALSE, xlab="",ylab="")
axis(2, 1:ncol(m), c("GM12878","K562","HepG2","H1"), las=2)
axis(1, 1:nrow(m), rownames(m), las=1)
for (y in 1:ncol(m))
	for (x in 1:nrow(m))
		text(x, y, format(m[x,y]*100,digits=2))

dev.off()


pdf("tissue_of_origin_20160817.heatmap.4tissues_encode.infer_best_subset.observed.pdf", paper="special", height=5, width=20)
par(oma=c(0, 0, 0, 0),mar=c(5, 5, 3, 1))
m<-t(apply(obs_ratio.no_overlap, 1, function(x)abs(x)/sum(abs(x),na.rm=T)))
image(1:nrow(m), 1:ncol(m), m, col = colorRampPalette(c("white","red"))(100), axes = FALSE, xlab="",ylab="")
axis(2, 1:ncol(m), c("GM12878","K562","HepG2","H1"), las=2)
axis(1, 1:nrow(m), rownames(m), las=1)
for (y in 1:ncol(m))
	for (x in 1:nrow(m))
		text(x, y, format(m[x,y]*100,digits=2))

dev.off()

pdf("tissue_of_origin_20160817.heatmap.4tissues_encode.infer_best_subset.RMSE.pdf", paper="special", height=5, width=10)
par(oma=c(0, 0, 0, 0),mar=c(4, 4, 3, 1))
barplot(rmses.no_overlap*100,ylab="RMSE on percentage",xlab="Samples",main="",ylim=c(0,20),names.arg=1:nrow(m))
abline(h=1,lty=2)
dev.off()

