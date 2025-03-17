library(DESeq2)

setwd("ColoRectal_Cancer/mNET-Seq/Replicates/Merged_ReadCounts/")
gdzie<-getwd()

## Load the readcount files
sampleFiles <- grep("RM-STARReadsPerGene.out.tab",list.files(gdzie),value=TRUE)

## Provide the names, condition, replicate number and order of the files to create sampleTable
sampleNames<-c( "HCEC_1CT-T4ph", "HCEC_1CT-Total", "HCT116-T4ph", "HCT116-Total","SW480-T4ph", "SW480-Total", "SW620-T4ph", "SW620-Total")
sampleCondition<-c("T4ph", "Total", "T4ph", "Total", "T4ph", "Total", "T4ph", "Total")
replicate<-c(1,1,1,1,1,1,1,1)
order<-c(1,2,3,4,5,6,7,8)

## SampleTable that has the information about the samples
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = factor(sampleCondition),replicate=factor(replicate), name=sampleNames)

##Create a DDSeq object
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = gdzie, design = ~ condition)
readCounts<-counts(dds, normalized = FALSE, replaced = FALSE)
readCountsSum<-colSums(readCounts)/1e6
dds2<-estimateSizeFactors(dds)
SF<-sizeFactors(dds2)
readCountsSumNorm<-readCountsSum/SF
readCountsSumNorm

### Load the genome coverage files(bedgraphs) for the normalisation
setwd("ColoRectal_Cancer/mNET-Seq/Replicates/Stranded_bedGraphs/")
gdzie<-getwd()
sampleFiles2 <- grep(".bedgraph.gz",list.files(gdzie),value=TRUE)

## Create another SampleTable with the bedGraphs
sampleNames2<-c("HCEC_1CT-T4ph-fwd", "HCEC_1CT-T4ph-rev", "HCEC_1CT-Total-fwd", "HCEC_1CT-Total-rev", "HCT116-T4ph-fwd", "HCT116-T4ph-rev", "HCT116-Total-fwd", "HCT116-Total-rev", "SW480-T4ph-fwd", "SW480-T4ph-rev", "SW480-Total-fwd", "SW480-Total-rev", "SW620-T4ph-fwd", "SW620-T4ph-rev", "SW620-Total-fwd", "SW620-Total-rev")
sampleCondition2<-c("T4ph", "T4ph", "Total", "Total", "T4ph", "T4ph", "Total", "Total", "T4ph", "T4ph", "Total", "Total", "T4ph", "T4ph", "Total", "Total")
replicate2<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
order2 <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
strand <-c("Fwd","Rev","Fwd","Rev","Fwd","Rev","Fwd","Rev","Fwd","Rev","Fwd","Rev","Fwd","Rev","Fwd","Rev")
SampleFactor<-rep(readCountsSumNorm, each = 2)
sampleTable2 <- data.frame(sampleName = sampleFiles2, fileName = sampleFiles2, condition = factor(sampleCondition2), replicate=factor(replicate2), name=sampleNames2, SampleFactor, strand)

chr<-c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM")

## bedtoolsBdgNorm is a function that normalise the genome coverage using the size-factor calculated above
bedtoolsBdgNorm<-function(x){
  f<-sampleTable2$SampleFactor[x]
  a<-read.table(sampleTable2[x,1], header=FALSE, stringsAsFactors=FALSE)
  a[,2]<-format(a[,2], scientific=FALSE)
  a[,3]<-format(a[,3], scientific=FALSE)
  if (sampleTable2$strand[x]=="Fwd") {
  a[,4]<-format(a[,4]/f,scientific=FALSE)} else { 
  a[,4]<-format(-a[,4]/f, scientific=FALSE)}
  a<-a[which(a[,1]%in% chr),]
  write.table(a ,file=paste(sampleTable2$name[x], "_norm.bedgraph", sep=""),quote=FALSE,sep="\t",row.names=FALSE, col.names=FALSE)
} 

for (i in 1:16) bedtoolsBdgNorm(i)
