library(DEXSeq)
library(GenomicRanges)
library(dplyr)
library(tidyr)

setwd("3Prime-Seq/IPF/polyAsites/PAS/PA_counts/CellLineWise_RObjects/1CTvsHCT116/")

### Load and process data (loops through fwd and rev files and assigns strand information and returns fwd/rev_data as lists because of multiple samples)
load_and_process_data <- function(file_pattern) {
  files <- list.files(pattern = file_pattern)
  lapply(files, function(f) {
    load(f)
    if (grepl("Fwd", f)) strand(binsBdgF) <- "+" else strand(binsBdgR) <- "-"
    if (grepl("Fwd", f)) binsBdgF else binsBdgR
  })
}

fwd_data <- load_and_process_data("FwdApa.RData")
rev_data <- load_and_process_data("RevApa.RData")

### Combines Forward and Reverse Lists => Combines the corresponding GRanges objects from Flist and Rlist into a single list
combined_data <- mapply(c, fwd_data, rev_data, SIMPLIFY = FALSE)

### Sort the combined data => Sorts the GRanges objects by sequence levels and genomic coordinates
combined_data_sorted <- lapply(combined_data, function(gr) {
    gr <- sortSeqlevels(gr)
    sort(gr)
})

### load the genes of interest (Major PAS annotation created as described in the methods of paper)
genes <- read.table("3Prime-Seq/APA_Analysis/1CT-major_PAS-SNR-6KbExt-3UTRExtended_wIds-AllPAS_Stringent_PCGs-Active.bed", header = TRUE, sep = "\t")
### Convert genes.bed file to Granges object for further downstream analysis 
gr_genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

assignGeneNames <- function(gr, genes) {
  # Find overlaps between 'gr' and 'genes'
  overlaps <- findOverlaps(gr, genes)
  
  # Initialize gene_id and gene_symbol columns as NA
  gr$gene_id <- rep(NA, length(gr))
  gr$gene_symbol <- rep(NA, length(gr))
  
  # Assign both gene_id and gene_symbol based on overlaps
  gr$gene_id[queryHits(overlaps)] <- genes$gene_id[subjectHits(overlaps)]
  gr$gene_symbol[queryHits(overlaps)] <- genes$gene_symbol[subjectHits(overlaps)]  # Populate gene_symbol
  
  # Create a data frame and filter out entries with missing gene_id
  df <- as.data.frame(gr) %>% filter(!is.na(gene_id)) %>% group_by(gene_id) %>% mutate(exon_number = paste0(gene_id, ":Exon", sprintf("%03d", row_number())), featureID = paste0("PAS", row_number())) %>% ungroup()
  
  # Convert back to GRanges object, keeping extra columns including gene_symbol
  result_gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  
  return(result_gr)
}

### Assigning geneNames to  combined data and creation of DEXSeq dataset and differential exon usage analysis
analyze_data <- function(data, genes) {
  annotated_data <- lapply(data, assignGeneNames, genes = genes)
  sum_matrix <- sapply(annotated_data, function(gr) gr$sum)

  featureIDs <- annotated_data[[1]]$featureID
  gene_ids <- annotated_data[[1]]$gene_id 
  df <- data.frame(sample = as.factor(c("1CT_1", "1CT_2", "1CT_3", "HCT116_2", "HCT116_3", "HCT116_4")), condition = as.factor(c("Control", "Control", "Control", "Tumor", "Tumor", "Tumor")), biolrep = as.factor(c("1", "2", "3", "1", "2", "3")))

  dataset <- DEXSeqDataSet(countData = sum_matrix, sampleData = df, design = ~ sample + exon + condition:exon, featureID = featureIDs, groupID = gene_ids)
  
  dataset <- dataset %>% estimateSizeFactors() %>% estimateDispersions() %>% testForDEU() %>% estimateExonFoldChanges(fitExpToVar="condition")
  
  # Save dispersion plot
  png("3Prime-Seq/APA_Analysis/Major_PAS-Annotation/1CTvsHCT116-Dispersion_PCGs-MP1A.png", width=8, height=6, units="in", res=200)
  plotDispEsts(dataset)
  dev.off()
  
  dxr <- DEXSeqResults(dataset)
  
  # Save MA plot
  png("3Prime-Seq/APA_Analysis/Major_PAS-Annotation/1CTvsHCT116-MA_Plot_PCGs-MP1A.png", width=8, height=6, units="in", res=200)
  plotMA(dxr, cex=0.8)
  dev.off()

  # Return both the DEXSeq results and the annotated data
  list(dxr_results = as.data.frame(dxr), annotated_data = annotated_data)
}

#results <- analyze_data(combined_data_sorted, genes_exp)
results <- analyze_data(combined_data_sorted, gr_genes)
### Access results
dxr <- results$dxr_results

# Check the structure of columns
sapply(dxr, is.list)

# Extract certain range of columns as one of the columns(genomicrange) in the dataframe is a list and that doesn't work with write.table function
colnames(dxr)
dxr1 <- cbind(dxr[,1:10],dxr[,12:17])

annotated_data <- results$annotated_data
gr_ctr <- annotated_data[[1]]
gr_ctr$sum <- NULL

gr_ctr1 <- as.data.frame(gr_ctr)
dxres <- cbind(dxr1, gr_ctr1[, c(1:5, 7)])

### Check few stats
length(unique(dxr$groupID))
table ( tapply( dxr$padj < 0.05, dxr$groupID, any))
table ( tapply( dxr$padj < 0.05, dxr$groupID, function(x) any(x, na.rm = TRUE)))

### Save Differential Exon usage results
write.table(dxres, "3Prime-Seq/APA_Analysis/Major_PAS-Annotation/1CTvsHCT116-DEXSeq_results_PCGs-MP1A.txt", sep="\t", quote=FALSE, row.names=FALSE)
dfr1 <- as.data.frame(dxr1)

### Classify the differential exon state as "shift" when a PAS padj < 0.05 orelse "no_shift" 
#group_classification <- tapply(dfr1$padj < 0.05, dfr1$groupID, function(x) if (any(x, na.rm = TRUE)) "shift" else "no_shift")

# Modified classification function that checks for both padj and log2FoldChange
group_classification <- tapply(dfr1$padj < 0.05 & !is.na(dfr1$log2fold_Tumor_Control), dfr1$groupID, function(x) if (any(x, na.rm = TRUE)) "shift" else "no_shift")
dim(group_classification)
dfr1$state <- group_classification[as.character(dfr1$groupID)]

# How many genes have one PAS?
df_temp <- dfr1 %>% group_by(groupID) %>% filter(n() ==1)
nrow(df_temp)

# Filter out genes that have more than one PAS
df2 <- dfr1 %>% group_by(groupID) %>% filter(n() > 1) %>% ungroup() %>% as.data.frame()
temp <- df2 %>% group_by(groupID, state) %>% summarize(n())
nrow(temp)

temp %>% group_by(state) %>% summarize(n())

n_distinct(df2[df2$state == "shift",]$groupID)
n_distinct(df2[df2$state == "no_shift",]$groupID)

### Save the temp in a file
write.table(temp, "3Prime-Seq/APA_Analysis/Major_PAS-Annotation/1CTvsHCT116-diff-APA_Shift_status_PCGs-MP1A.txt", row.names=FALSE, quote=FALSE, sep="\t")

### APA Analysis
# First part of analysis: Extract the top 2 PAS and assign its localization based on strand and start as mentioned in MOl Cell Paper
apa_analysis_part1 <- function(df, gr) {
  df %>% filter(state == "shift") %>% mutate(count_sum = rowSums(select(., starts_with("countData.")))) %>%
    group_by(groupID) %>% arrange(padj, desc(count_sum), .by_group = TRUE) %>% slice_head(n = 2) %>%
    left_join(as.data.frame(gr_ctr), by = c("groupID" = "gene_id", "featureID")) %>% filter(n_distinct(strand) == 1) %>%
    mutate(localization = case_when(
      strand == "+" & start == min(start) ~ "proximal",
      strand == "+" & start != min(start) ~ "distal",
      strand == "-" & start == max(start) ~ "proximal",
      strand == "-" & start != max(start) ~ "distal"
    )) %>% ungroup()
}

# Second part of analysis: calculate ratios and determine direction
apa_analysis_part2 <- function(data) {
  data %>% group_by(groupID) %>% summarize(
      ratio_Tumor = Tumor[localization == "distal"] / Tumor[localization == "proximal"],
      ratio_Control = Control[localization == "distal"] / Control[localization == "proximal"]
    ) %>% mutate(direction = case_when(
      ratio_Tumor > ratio_Control ~ "distal",
      ratio_Tumor < ratio_Control ~ "proximal",
      ratio_Tumor == ratio_Control ~ "DP"
    ))
}

intermediate_results <- apa_analysis_part1(df2, gr_ctr)

write.table(intermediate_results, "3Prime-Seq/APA_Analysis/Major_PAS-Annotation/1CTvsHCT116-diff-APA-results_PCGs-MP1A.txt", sep="\t", quote=FALSE, row.names=FALSE)
  
final_results <- apa_analysis_part2(intermediate_results)
write.table(final_results, "3Prime-Seq/APA_Analysis/Major_PAS-Annotation/1CTvsHCT116-diff-APA_Ratio_PCGs-MP1A.txt", sep="\t", quote=FALSE, row.names=FALSE)

final_results %>% group_by(direction) %>% summarize(n = n())
