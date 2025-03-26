library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("3Prime-Seq/IPF/polyAsites/PAS/Raw_PAS/")

FwdCovGr<-list.files(pattern="bam_raw_PAS_Fwd.RData")  ## lists only PAS, not control files
RevCovGr<-list.files(pattern="bam_raw_PAS_Rev.RData")
Flist<-list(1:length(FwdCovGr))
Rlist<-list(1:length(RevCovGr))

for (i in 1:length(FwdCovGr)) {
	load(FwdCovGr[i])
	Flist[i]<-grr
}

FwdCovTot<-do.call("c", Flist) # merges all list elements into on GRanges object
f<-length(FwdCovTot)/1e6
FwdCovTotS<-sort(FwdCovTot)
FwdCovTotScov<-coverage(FwdCovTotS)
FwdCovTotGfin<-GRanges(FwdCovTotScov)

covG<-FwdCovTotGfin
save(covG,file="all_PAS_covG_Fwd.RData")
end(covG)<-end(covG)+1
covG<-trim(covG)
a<-as.data.frame(covG)
a<-a[,c(1:3,6)]
a[,2]<-format(a[,2], scientific=FALSE)
a[,3]<-format(a[,3], scientific=FALSE) ## 1 added for proper display in UCSC
a[,4]<-format(a[,4]/f, scientific=FALSE)
write.table(a ,file="all_PAS_covG_Fwd.bedGraph",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)


for (i in 1:length(RevCovGr)) {
	load(RevCovGr[i])
	Rlist[i]<-grr
}

RevCovTot<-do.call("c", Rlist) # merges all list elements into on GRanges object
f<-length(RevCovTot)/1e6
RevCovTotS<-sort(RevCovTot)
RevCovTotScov<-coverage(RevCovTotS)
RevCovTotGfin<-GRanges(RevCovTotScov)

covG<-RevCovTotGfin
save(covG,file="all_PAS_covG_Rev.RData")
end(covG)<-end(covG)+1
covG<-trim(covG)
a<-as.data.frame(covG)
a<-a[,c(1:3,6)]
a[,2]<-format(a[,2], scientific=FALSE)
a[,3]<-format(a[,3], scientific=FALSE) ## 1 addef for proper display in UCSC
a[,4]<-format(a[,4]/f, scientific=FALSE)
write.table(a ,file="all_PAS_covG_Rev.bedGraph",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)


### preparing set of all PAS for samples from both events

load("all_PAS_covG_Fwd.RData")
FwdPile<-covG
FwdPile30<-FwdPile[which(FwdPile$score>30)]
FwdPile30r<-reduce(FwdPile30)
table(width(FwdPile30r))

load("all_PAS_covG_Rev.RData")
RevPile<-covG
RevPile30<-RevPile[which(RevPile$score>30)]
RevPile30r<-reduce(RevPile30)
table(width(RevPile30r))

### finding local maxima round 1
# find all coverage intervals with signal >30 (cov from 12 samples)
# find maximal PAS signal within local interval & make 30nt window around it 

FwdGlist<-list(1:length(FwdPile30r))
for (i in 1:length(FwdPile30r)) {
	x<-subsetByOverlaps(FwdPile30,FwdPile30r[i])
	x_max<-x[which(x$score==max(x$score))][1]
	x_max<-resize(x_max,30,fix="center")
	FwdGlist[i]<-x_max
}
save(FwdGlist,file="PAS_max1_Fwd.RData")

RevGlist<-list(1:length(RevPile30r))
for (i in 1:length(RevPile30r)) {
	x<-subsetByOverlaps(RevPile30,RevPile30r[i])
	x_max<-x[which(x$score==max(x$score))][1]
	x_max<-resize(x_max,30,fix="center")
	RevGlist[i]<-x_max
}
save(RevGlist,file="PAS_max1_Rev.RData")

### Find local maxima round 2
# Create reduced GRanges object from list of 30nt windows centered at the local max PAS signal
# Overlap its intervals with all coverage intervals with signal >30 (cov from 36 samples)
# Find maximum & make 30 nt window around it

FwdG1<-do.call("c", FwdGlist)
FwdG1r<-reduce(FwdG1)

RevG1<-do.call("c", RevGlist)
RevG1r<-reduce(RevG1)

FwdGlist2<-list(1:length(FwdG1r))
for (i in 1:length(FwdG1r)) {
	x<-subsetByOverlaps(FwdPile30,FwdG1r[i])
	x_max<-x[which(x$score==max(x$score))][1]
	x_max<-resize(x_max,30,fix="center")
	FwdGlist2[i]<-x_max
}
save(FwdGlist2,file="PAS_max2_Fwd.RData")

RevGlist2<-list(1:length(RevG1r))
for (i in 1:length(RevG1r)) {
	x<-subsetByOverlaps(RevPile30,RevG1r[i])
	x_max<-x[which(x$score==max(x$score))][1]
	x_max<-resize(x_max,30,fix="center")
	RevGlist2[i]<-x_max
}
save(RevGlist2,file="PAS_max2_Rev.RData")

FwdG2<-do.call("c", FwdGlist2)
FwdG2Df<-as.data.frame(FwdG2)
write.table(FwdG2Df[,1:3] ,file="PAS_local_max2_Fwd.bed" ,quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)

RevG2<-do.call("c", RevGlist2)
RevG2Df<-as.data.frame(RevG2)
write.table(RevG2Df[,1:3] ,file="PAS_local_max2_Rev.bed" ,quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
