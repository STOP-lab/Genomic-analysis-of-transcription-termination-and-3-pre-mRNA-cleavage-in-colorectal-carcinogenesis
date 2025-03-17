library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

## Listing stranded Fwd/Rev BAM files with alignments
Fwd <- list.files(path = "Replicates/Stranded_BAMs/", pattern = "fwd", full.names = TRUE)
Rev <- list.files(path = "Replicates/Stranded_BAMs/", pattern = "rev", full.names = TRUE)

## Extractor is a function that is used to extract the last transcribed nucleotide, i.e., the last nucleotide of a read and it also includes two sub-functions Bam2bedgraph_Fwd, Bam2bedgraph_Fwd that converts the BAM file to bedgraph with the last nucleotide extracted.

extractor <- function(i) {
    tryCatch({
        # Forward strand
        F <- readGAlignments(Fwd[i])
        F <- keepStandardChromosomes(F, pruning.mode = "coarse")
        F <- dropSeqlevels(F, "chrM", pruning.mode = "coarse")
        Fg <- GRanges(F)
        start(Fg) <- end(Fg)  # Last nucleotide of the forward reads
        
        # Reverse strand
        R <- readGAlignments(Rev[i])
        R <- keepStandardChromosomes(R, pruning.mode = "coarse")
        R <- dropSeqlevels(R, "chrM", pruning.mode = "coarse")
        Rg <- GRanges(R)
        end(Rg) <- start(Rg)  # Last nucleotide of the reverse reads
        
        # Save forward and reverse data
        Bam2bedgraph_Fwd(Fg, basename(Fwd[i]))
        Bam2bedgraph_Rev(Rg, basename(Rev[i]))
    }, error = function(e) {
        message(paste("Error processing file index", i, ":", e$message))
    })
}

Bam2bedgraph_Fwd <- function(gr, name) {
    Path <- "Replicates/Stranded_bedGraphs/SNR/"
    grr <- gr
    start(grr) <- end(gr) - 1
    
    covG <- GRanges(coverage(grr), seqinfo = seqinfo(Hsapiens))
    save(covG, file = file.path(Path, paste0(name, "-SNR_raw_Fwd.RData")))
    
    a <- as.data.frame(covG)
    a <- a[, c(1:3, 6)]
    a[, 2] <- format(a[, 2], scientific = FALSE)
    a[, 3] <- format(a[, 3], scientific = FALSE)
    write.table(a, file = file.path(Path, paste0(name, "-SNR_raw_Fwd.bedGraph")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

Bam2bedgraph_Rev <- function(gr, name) {
    Path <- "Replicates/Stranded_bedGraphs/SNR/"
    grr <- gr
    end(grr) <- start(gr) + 1
    
    covG <- GRanges(coverage(grr), seqinfo = seqinfo(Hsapiens))
    save(covG, file = file.path(Path, paste0(name, "-SNR_raw_Rev.RData")))
    
    a <- as.data.frame(covG)
    a <- a[, c(1:3, 6)]
    a[, 2] <- format(a[, 2], scientific = FALSE)
    a[, 3] <- format(a[, 3], scientific = FALSE)
    write.table(a, file = file.path(Path, paste0(name, "-SNR_raw_Rev.bedGraph")), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

for (i in 1:length(Fwd)) {
    extractor(i)
}
