library(DESeq2)
setwd("ColoRectal_Cancer/POINT-Seq/Replicates/Merged_ReadCounts/")
gdzie<-getwd()

sampleFiles <- grep("STARReadsPerGene.out.tab",list.files(gdzie),value=TRUE)
sampleFiles

sampleNames<-c("SW480_DMSO", "SW480_JTE", "SW620_DMSO", "SW620_JTE")

sampleCondition<-c("Ctrl", "Treatment", "Ctrl", "Treatment")
replicate<-c(1,1,1,1)
order<-c(1,2,3,4)
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

setwd("ColoRectal_Cancer/POINT-Seq/Replicates/Stranded_bedGraphs/")
gdzie<-getwd()
sampleFiles2 <- grep(".bedgraph.gz",list.files(gdzie),value=TRUE)
sampleFiles2

sampleNames2 <-c("SW480_DMSO-fwd", "SW480_DMSO-rev", "SW480_JTE-fwd",  "SW480_JTE-rev", "SW620_DMSO-fwd", "SW620_DMSO-rev", "SW620_JTE-fwd", "SW620_JTE-rev")

sampleCondition2<-c("Ctrl", "Ctrl", "Treatment", "Treatment", "Ctrl", "Ctrl", "Treatment", "Treatment")

replicate2<-c(1,1,1,1,1,1,1,1)
order2 <-c(1,2,3,4,5,6,7,8)
strand <-c("Fwd","Rev","Fwd","Rev","Fwd","Rev","Fwd","Rev")

SampleFactor<-rep(readCountsSumNorm, each = 2)
SampleFactor

sampleTable2 <- data.frame(sampleName = sampleFiles2, fileName = sampleFiles2, condition = factor(sampleCondition2), replicate=factor(replicate2), name=sampleNames2, SampleFactor, strand)
sampleTable2

chr<-c(paste("chr",1:22,sep=""),"chrX", "chrY", "chrM")

bedtoolsBdgNorm<-function(x){
  Path <- "ColoRectal_Cancer/POINT-Seq/Replicates/Normalised_bedGraphs/"
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

for (i in 1:8) bedtoolsBdgNorm(i)

