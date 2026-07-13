library(DESeq2)

setwd("/dysk2/groupFolders/deepshika/mNET-Seq/Pancreatic_Cells/Replicates/Merged_ReadCounts/")
gdzie<-getwd()

sampleFiles <- grep("STAR-SEReadsPerGene.out.tab",list.files(gdzie),value=TRUE)
sampleFiles

sampleNames<-c("BxPC3-T4ph-SNR", "MiaPaCa2-T4ph-SNR", "Panc1-T4ph-SNR")

sampleCondition<-c("WT", "Mutated", "Mutated")
replicate<-c(1,1,1)
order<-c(1,2,3)
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = factor(sampleCondition),replicate=factor(replicate), name=sampleNames)
sampleTable

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = gdzie, design = ~ condition)
readCounts<-counts(dds, normalized = FALSE, replaced = FALSE)
readCountsSum<-colSums(readCounts)/1e6
readCountsSum

dds2<-estimateSizeFactors(dds)
SF<-sizeFactors(dds2)
SF

readCountsSumNorm<-readCountsSum/SF
readCountsSumNorm

### Single Nucleotide Resolution files processing
#setwd("/dysk2/groupFolders/deepshika/mNET-Seq/Pancreatic_Cells/Replicates/Merged_BAMs/SNR/Stranded_bedGraphs/")
setwd("/dysk2/groupFolders/deepshika/mNET-Seq/Pancreatic_Cells/Replicates/Merged_BAMs/SNR/Stranded_bedGraphs/Markdup/")

gdzie<-getwd()
sampleFiles2 <- grep(".bedgraph",list.files(gdzie),value=TRUE)
sampleFiles2

#sampleNames2<-c("BxPC3-T4ph-dd-SNR-fwd", "BxPC3-T4ph-dd-SNR-rev", "MiaPaCa2-T4ph-dd-SNR-fwd", "MiaPaCa2-T4ph-dd-SNR-rev", "Panc1-T4ph-dd-SNR-fwd", "Panc1-T4ph-dd-SNR-rev")
sampleNames2<-c("BxPC3-T4ph-SNR-fwd", "BxPC3-T4ph-SNR-rev", "MiaPaCa2-T4ph-SNR-fwd", "MiaPaCa2-T4ph-SNR-rev", "Panc1-T4ph-SNR-fwd", "Panc1-T4ph-SNR-rev")

sampleCondition2<-c("WT", "WT", "Mutated", "Mutated", "Mutated", "Mutated")
replicate2<-c(1,1,1,1,1,1)
order2 <-c(1,2,3,4,5,6)
strand <-c("Fwd","Rev","Fwd","Rev","Fwd","Rev")

SampleFactor<-rep(readCountsSumNorm, each = 2)
SampleFactor

sampleTable2 <- data.frame(sampleName = sampleFiles2, fileName = sampleFiles2, condition = factor(sampleCondition2), replicate=factor(replicate2), name=sampleNames2, SampleFactor, strand)
sampleTable2

chr<-c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM")

bedtoolsBdgNorm<-function(x){
#  Path <- "/dysk2/groupFolders/deepshika/mNET-Seq/Pancreatic_Cells/Replicates/Merged_BAMs/SNR/Normalised_bedGraphs/"
  Path <- "/dysk2/groupFolders/deepshika/mNET-Seq/Pancreatic_Cells/Replicates/Merged_BAMs/SNR/Normalised_bedGraphs/Markdup/"
  f <- sampleTable2$SampleFactor[x]
  a <- read.table(sampleTable2[x, 1], header = FALSE, stringsAsFactors = FALSE)
  a[, 2] <- as.integer(a[, 2]) 
  a[, 3] <- as.integer(a[, 3])
  if (sampleTable2$strand[x] == "Fwd") {
   a[, 4] <- a[, 4] / f} else { 
   a[, 4] <- -a[, 4] / f}
  # Round to desired decimals
  a[, 4] <- round(a[, 4], digits = 6)   
  a <- a[which(a[, 1] %in% chr), ]
  write.table(a, file = file.path(Path, paste(sampleTable2$name[x], "_norm.bedgraph", sep = "")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

for (i in 1:6) bedtoolsBdgNorm(i)
