library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

Fwd_apa<-read.table("3Prime-Seq/IPF/polyAsites/PAS/APA_Resources/PAS_local_max2_Fwd.bed", header=FALSE, stringsAsFactors=FALSE)
Rev_apa<-read.table("3Prime-Seq/IPF/polyAsites/PAS/APA_Resources/PAS_local_max2_Rev.bed", header=FALSE, stringsAsFactors=FALSE)

Fwd_apaG<-GRanges(Fwd_apa[,1],IRanges(Fwd_apa[,2], Fwd_apa[,3]))
Rev_apaG<-GRanges(Rev_apa[,1],IRanges(Rev_apa[,2], Rev_apa[,3]))

setwd("3Prime-Seq/IPF/polyAsites/PAS/Raw_PAS/")

FwdPAgr<-list.files(pattern="bam_raw_PAS_Fwd.RData")  ## lists only PAs, not control files
RevPAgr<-list.files(pattern="bam_raw_PAS_Rev.RData")

ApaQuant <- function(i) {
    # Print the current index and file paths for debugging
    cat("Processing index:", i, "\n")
    cat("Loading forward PA file:", FwdPAgr[i], "\n")
    cat("Loading reverse PA file:", RevPAgr[i], "\n")
    
    # Check if the file paths are not NA
    if (is.na(FwdPAgr[i]) || is.na(RevPAgr[i])) {
        cat("Skipping index", i, "due to NA file path\n")
        return(NULL)
    }

    # Load forward and reverse PA events
    load(FwdPAgr[i])
    F <- grr
    load(RevPAgr[i])
    R <- grr

    # Calculate coverage for forward and reverse PA events
    covF <- coverage(F)
    covR <- coverage(R)

    # Function to calculate the sum of coverage within bins
    binnedAverage <- function(bins, numvar, mcolname) {
        stopifnot(is(bins, "GRanges"))
        stopifnot(is(numvar, "RleList"))
        stopifnot(identical(seqlevels(bins), names(numvar)))
        bins_per_chrom <- split(ranges(bins), seqnames(bins))
        means_list <- lapply(names(numvar), function(seqname) {
            views <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
            viewSums(views)
        })
        new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
        mcols(bins)[[mcolname]] <- new_mcol
        bins
    }

    # Apply binnedAverage to forward and reverse APA sites
    binsBdgF <- binnedAverage(Fwd_apaG, covF[1:24], "sum")
    binsBdgR <- binnedAverage(Rev_apaG, covR[1:24], "sum")

    # Save the binned average results
    save(binsBdgF, file=paste("/dysk2/groupFolders/deepshika/ColoRectal_Cancer/3Prime-Seq/IPF/polyAsites/PAS/PA_counts/", FwdPAgr[i], "_FwdApa.RData", sep=""))
    save(binsBdgR, file=paste("/dysk2/groupFolders/deepshika/ColoRectal_Cancer/3Prime-Seq/IPF/polyAsites/PAS/PA_counts/", RevPAgr[i], "_RevApa.RData", sep=""))

    # Convert results to data frames and adjust start/end positions
    binsBdgFdf <- as.data.frame(binsBdgF)
    binsBdgFdf[,2] <- format(binsBdgFdf[,2] + 14, scientific=FALSE)
    binsBdgFdf[,3] <- format(binsBdgFdf[,3] - 14, scientific=FALSE)
    binsBdgFdf[,7] <- binsBdgFdf[,6] / (length(F) / 1e6) # Normalize to RPM

    binsBdgRdf <- as.data.frame(binsBdgR)
    binsBdgRdf[,2] <- format(binsBdgRdf[,2] + 14, scientific=FALSE)
    binsBdgRdf[,3] <- format(binsBdgRdf[,3] - 14, scientific=FALSE)
    binsBdgRdf[,7] <- -binsBdgRdf[,6] / (length(R) / 1e6) # Normalize to RPM and negate

    # Write results to bedGraph files
    write.table(binsBdgFdf[,c(1:3,7)], file=paste("3Prime-Seq/IPF/polyAsites/PAS/PA_counts/", FwdPAgr[i], "_FwdApaCounts_rpmC.bedGraph", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(binsBdgRdf[,c(1:3,7)], file=paste("3Prime-Seq/IPF/polyAsites/PAS/PA_counts/", RevPAgr[i], "_RevApaCounts_rpmC.bedGraph", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

# Run the function for each index within the valid range
for (i in 1:length(FwdPAgr)) {
    ApaQuant(i)
}
