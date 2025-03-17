#### filtering of Lexogen QuantSeq Rev data to keep only genuine polyadenylation sites & remove signal internal primed reads 

### Filtering out of internal priming events; modified from Derti et al. Genome Res. 2012. 22: 1173-1183 (consider upadting of APAs from https://academic.oup.com/nar/article/48/D1/D174/5588346)  

# Overview:
# - find reads, where the genomic sequence directly downstream contains a stretch of As: potential internal priming (Ts on the - strand)               
#     -- > 6/10 bases are A
#     --6 consecutive bases in the 10nt window are A
# 
# 2. exclude from the mask:
# - same potential sites that overlap +-20bp a 3'end from GENCODE annptaion (gtf file): 
# - exclude polyA sites from Derti et al. 2012 (downloaded from UCSC table browser, available also form GEO)
# [optional 2: looking for overlap with known 3'ends might be good enough, otherwise exclude also potential internal priming sites containing PAS signal (see below) -40 to -10 from 3' end]

### 1. creating new GTF files with annotation for beginnings and end of the whole genes
### 2. merging all APA events form different tissues from Derti et al. 2012; data in hg19 neds to be lifted to hg38

### merging all APA events form different tissues from Derti et al. 2012; data in hg19 neds to be lifted to hg38

R

library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library("rtracklayer")
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/polyAsites/")

files<- list.files()

Fwd<-list.files(pattern="_F")
Rev<-list.files(pattern="_R")

Flist<-list(1:9)
Rlist<-list(1:9)

for (i in 1:9) {
	p<-read.table(Fwd[i], header=FALSE, skip=1, stringsAsFactors=FALSE)
	pG<-GRanges(p$V1,IRanges(p$V2,p$V3))
	pG<-keepStandardChromosomes(pG, species="Homo sapiens", pruning.mode="coarse")
	pG<-resize(pG, 20, fix="center", use.names=TRUE, ignore.strand=TRUE)  ### change of window size form 40 to 20
	pG<-sortSeqlevels(pG)
	pG<-sort(pG)
	Flist[i]<-pG
} 

FwdTSE<-do.call("c", Flist)
FwdTSEr<-reduce(FwdTSE)
export(FwdTSEr, "polyadenylationSites_merged_Derti_GenRes_2012_FwdTSEr_7jan24.bed")

for (i in 1:9) {
	p<-read.table(Rev[i], header=FALSE, skip=1, stringsAsFactors=FALSE)
	pG<-GRanges(p$V1,IRanges(p$V2,p$V3))
	pG<-keepStandardChromosomes(pG, species="Homo sapiens", pruning.mode="coarse")
	pG<-resize(pG, 20, fix="center", use.names=TRUE, ignore.strand=TRUE) ### change of window size form 40 to 20
	pG<-sortSeqlevels(pG)
	pG<-sort(pG)
	Rlist[i]<-pG
} 

RwdTSE<-do.call("c", Rlist)
RwdTSEr<-reduce(RwdTSE)
export(RwdTSEr, "polyadenylationSites_merged_Derti_GenRes_2012_RwdTSEr_7jan24.bed")

save(FwdTSEr, file="FwdTSEr.RData")
save(RwdTSEr, file="RwdTSEr.RData")

# liftOver hg19 hg38
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38) 

hg19ToHg38<-import.chain("/home/micgdu/Analysis/utilities/UCSC/chainFiles/human/hg19ToHg38.over.chain")

LiftGR<-function(xg, chain) {
	# GRanges object as input
	cur19 = liftOver(xg, chain)
	cur20<-unlist(cur19, recursive = TRUE, use.names = TRUE)
	cur21<-as.data.frame(cur20)
	cur21[,2]<-format(cur21[,2],scientififc=FALSE)
	cur21[,3]<-format(cur21[,3],scientififc=FALSE)
	return(cur21[,c(1:3)])
}

FwdTSEr_hg38<-LiftGR(FwdTSEr, hg19ToHg38)
RwdTSEr_hg38<-LiftGR(RwdTSEr, hg19ToHg38)

write.table(FwdTSEr_hg38, file="polyadenylationSites_merged_hg38_Derti_GenRes_2012_FwdTSEr_7jan24.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names = TRUE)
write.table(RwdTSEr_hg38, file="polyadenylationSites_merged_hg38_Derti_GenRes_2012_RwdTSEr_7jan24.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names = TRUE)
# bed files uploaded to the UCSC with the header:
# track name="T_mask4c" description="T_mask all TxA in a row & 10nt >6T which rae not 3'end or exp APAs " color=100,0,50 

save(FwdTSEr_hg38, file="FwdTSEr_hg38.RData")
save(RwdTSEr_hg38, file="RwdTSEr_hg38.RData")

### genome patternes and Gencode 3' ends

library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38) 

setwd("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter")
genes<-read.table("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/GTF_whole_genes_hg38_geneCodeV42_7jan24.gtf",header=FALSE,stringsAsFactors=FALSE)

# 1) generating 20bp widows around 3' ends
genesG<-GRanges(genes$V1,IRanges(genes$V4,genes$V5), strand=genes$V7)

genesEndPlus<-genesG[which(strand(genesG)=="+")]
start(genesEndPlus)<-end(genesEndPlus)-10
end(genesEndPlus)<-end(genesEndPlus)+10

genesEndMinus<-genesG[which(strand(genesG)=="-")]
end(genesEndMinus)<-start(genesEndMinus)+10
start(genesEndMinus)<-start(genesEndMinus)-10

genesEndPlusDf<-as.data.frame(genesEndPlus)
genesEndMinusDf<-as.data.frame(genesEndMinus)

write.table(genesEndPlusDf[,1:3] ,file="genesEndPlus10.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
write.table(genesEndMinusDf[,1:3] ,file="genesEndMinus10.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)

# 2) find "AAAAAA" stretch

hs<-getSeq(Hsapiens)
hs<-hs[1:25]

A6 <- BString("AAAAAA")   ### Fwd
T6 <- BString("TTTTTT")   ### Rev

#matchPattern(pattern, subject)
#matchPattern(pattern, subject)

A6match<-vmatchPattern(A6, hs)
A6matchR<-unlist(A6match, recursive=TRUE, use.names=TRUE)
A6matchG<-GRanges(names(A6matchR),A6matchR)
A6matchGr<-reduce(A6matchG)

T6match<-vmatchPattern(T6, hs)
T6matchR<-unlist(T6match, recursive=TRUE, use.names=TRUE)
T6matchG<-GRanges(names(T6matchR),T6matchR)
T6matchGr<-reduce(T6matchG)

# 3) find windows with >6 A/T in 10 bp windows

val<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20", "21","22","X","Y","M")
val<-paste("chr",val,sep="")

la<-list(1:25)
for (i in c(1:25)) {
	Aper10w<-letterFrequencyInSlidingView(hs[[i]],10, "A")
Aper10wInd<-cbind(1:length(Aper10w),Aper10w)
Aper10wIndLim<-Aper10wInd[which(Aper10wInd[,2]>6),]
Aper10wIndLimG<-GRanges(paste("chr",i,sep=""),IRanges(Aper10wIndLim[,1],Aper10wIndLim[,1]+9))
Aper10wIndLimGr<-reduce(Aper10wIndLimG)
la[i]<-Aper10wIndLimGr
}

A10<-do.call("c", la)
A10<-renameSeqlevels(A10, val)

lt<-list(1:25)
for (i in c(1:25)) {
	Aper10w<-letterFrequencyInSlidingView(hs[[i]],10, "T")
Aper10wInd<-cbind(1:length(Aper10w),Aper10w)
Aper10wIndLim<-Aper10wInd[which(Aper10wInd[,2]>6),]
Aper10wIndLimG<-GRanges(paste("chr",i,sep=""),IRanges(Aper10wIndLim[,1],Aper10wIndLim[,1]+9))
Aper10wIndLimGr<-reduce(Aper10wIndLimG)
lt[i]<-Aper10wIndLimGr
}

T10<-do.call("c", lt)
T10<-renameSeqlevels(T10, val)

# making final mask: Fwd: A6 + A10 - genesEndPlus
# 			 Rev  T6 + T10 - genesEndMinus

A_both<-reduce(c(A6matchGr,A10))
T_both<-reduce(c(T6matchGr,T10))

A_mask<-setdiff(A_both,genesEndPlus, ignore.strand=TRUE)
T_mask<-setdiff(T_both,genesEndMinus, ignore.strand=TRUE)

A_maskDf<-as.data.frame(A_mask)
T_maskDf<-as.data.frame(T_mask)

A_10Df<-as.data.frame(A10)
T_10Df<-as.data.frame(T10)

write.table(A_10Df[,1:3] ,file="A_10.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
write.table(T_10Df[,1:3] ,file="T_10.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)

write.table(A_maskDf[,1:3] ,file="A_mask.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
write.table(T_maskDf[,1:3] ,file="T_mask.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)

load("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/polyAsites/FwdTSEr_hg38.RData")
load("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/polyAsites/RwdTSEr_hg38.RData")

FwdTSEr<-FwdTSEr_hg38
RwdTSEr<-RwdTSEr_hg38

RwdTSErG<-makeGRangesFromDataFrame(RwdTSEr)
FwdTSErG<-makeGRangesFromDataFrame(FwdTSEr)

A_mask<-keepStandardChromosomes(A_mask, species="Homo_sapiens", pruning.mode="coarse")
T_mask<-keepStandardChromosomes(T_mask, species="Homo_sapiens", pruning.mode="coarse")

A_mask2<-setdiff(A_mask,FwdTSErG, ignore.strand=TRUE)
T_mask2<-setdiff(T_mask,RwdTSErG, ignore.strand=TRUE)

A_maskDf2<-as.data.frame(A_mask2)
T_maskDf2<-as.data.frame(T_mask2)

save(A_mask2,file="A_mask4.RData")
save(T_mask2,file="T_mask4.RData")

write.table(A_maskDf2[,1:3] ,file="A_mask_polyAseq_hg38.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
write.table(T_maskDf2[,1:3] ,file="T_mask_polyAseq_hg38.bed",quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)

#load("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/A_mask4.RData")
#load("/home/micgdu/Analysis/Kinga/3endQuantSec/IntPrFilter/T_mask4.RData")



