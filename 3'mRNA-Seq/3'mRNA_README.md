# 3'mRNA-Seq Analysis
- 3'mRNA-Seq has been performed on Colorectal cells. Please refer to the script for all sample processing [1.3Prime-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/1.3Prime-Seq_processing.py)
- [Samples_3'mRNA](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/Samples_3'mRNA) - contains the name of fastq files used to process all samples at once
# Folder structure for analysis
```
mkdir 3Prime-Seq/FastQ 3Prime-Seq/QC 3Prime-Seq/Trimmed 3Prime-Seq/Trimmed/PolyA-T/ 3Prime-Seq/Aligned mkdir 3Prime-Seq/Aligned/ReadCounts/
mkdir 3Prime-Seq/QC/FastQC 3Prime-Seq/QC/TrimmedQC 3Prime-Seq/QC/TrimmedQC/PolyA-T 3Prime-Seq/QC/Aligned
mkdir 3Prime-Seq/Replicates 3Prime-Seq/Aligned/Stranded_BAMs/ 3Prime-Seq/Aligned/Stranded_bedGraphs/
mkdir 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/ 3Prime-Seq/Aligned/Stranded_bigWigs/
mkdir 3Prime-Seq/IPF 3Prime-Seq/Aligned/IPF/PolyAsites 3Prime-Seq/IPF/PAS 3Prime-Seq/IPF/Raw_PAS
mkdir 3Prime-Seq/IPF/PAs_counts 3Prime-Seq/APA-Analysis
mkdir 
```  
# 1. DeMultiplexing and FastQ format conversion
- Raw .cbcl files to fastq conversion and demultiplexing  was done using custom-created SampleSheet
- --no-lane-splitting Make sure that the reads of one sample from different lanes don't produce multiple fastq files but would end up in one FastQ file
```  	
bcl2fastq -R 231030_A01680_0155_AHMN5MDRX3/ -o 231030_A01680_0155_AHMN5MDRX3/3Prime-Seq/FastQ/ --no-lane-splitting -p 26
```
# 2. Quality check
```
fastqc 3Prime-Seq/FastQ/*.gz -o 3Prime-Seq/QC/FastQC/
```
# 3. Adapter, low-quality reads, and PolyA/T_tails removal
 - [1.3Prime-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/1.3Prime-Seq_processing.py)
 - Move the trimmed reports to the quality check reports to main the folder structure
```	
mv Trimmed/*.zip QC/TrimmedQC/
mv Trimmed/*.html QC/TrimmedQC/
mv Trimmed/*.txt QC/TrimmedQC/
mv Trimmed/*.zip QC/TrimmedQC/PolyA-T
mv Trimmed/*.html QC/TrimmedQC/PolyA-T
mv Trimmed/*.txt QC/TrimmedQC/PolyA-T
```
# 4. Alignment 
 - STAR aligner, human genome (hg38)
 - Parameters for STAR alignment were adapted from lexogen (https://github.com/Lexogen-Tools/quantseqpool_analysis) along with the generation of read counts
 - [1.3Prime-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/1.3Prime-Seq_processing.py)
# 5. Generate strand-specific BAMs and genome coverage files(bedGrahps)
- Forward Strand reads - 5' to 3' direction; use the flag -f 16
- Reverse Strand reads - 3' to 5' direction; use the flag -F 16
  -f => includes, -F => excludes;  16 => reads mapping to reverse strand
- [1.3Prime-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/1.3Prime-Seq_processing.py)
# 6. Normalisation of strand-specific reads
- Estimates size factors in Deseq2 and normalisation factors (per million factors)
- Generates normalised bedgraphs
- [2.Normalised_bedGraphs.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/2.Normalised_bedGraphs.R)
- Move Normalised bedgraphs to their folder
```
mv 3Prime-Seq/Aligned/Stranded_bedGraphs/*.bedGraph 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/
mv 3Prime-Seq/Aligned/Stranded_bedGraphs/*.bedGraph 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/
```
- Sort normalised bedgraphs and generate bigWigs
```
for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do sort -k1,1 -k2,2n "$f" -o "$f"; done
for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do /home/micgdu/kentutils/bedGraphToBigWig "$f" hg38_chromsizes.genome "$f.bw"; done
for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph.bw/.bw/g)"; done
mv 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw 3Prime-Seq/Aligned/Stranded_bigWigs/
```
# 7. Filtering Internal Priming Events
- Uses two scripts [3.InternalPrimming_Mask.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/3.InternalPrimming_Mask.R); [4.InternalPrimming_Filtering.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/4.InternalPrimming_Filtering.R)
- [3.InternalPrimming_Mask.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/3.InternalPrimming_Mask.R) : Generates the genomic mask that includes consecutive 6A's/T's and or >6A/T in a 10 nucleotide window excluding the gene 3'ends and validates experimental polyadenylation sites
- [4.InternalPrimming_Filtering.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/4.InternalPriming_Filtering.R) - Script used to filter out the internal priming reads from the strand-specific bams based on the crude genomic mask generated

# 8. Polyadenylation site (PAS) Detection and Quantification
- [5.PAS_Calling.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/5.PAS_Calling.R)
  Generate density profiles: Strand-specific, depth-normalized 3’ mRNA-seq coverage profiles were created for each strand across the genome
- [6.PAS_Quantification.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/6.PAS_Quantification.R)
- Identifying PAS: A sliding 30-nucleotide window summed signal values, selecting windows above a threshold (30)
                   The strongest peak within each window was identified, and overlapping regions were removed to refine PAS detection
                   The detected PAS were used to extract sample-specific polyadenylation sites
  
- Rename and relocate files after PAS Quantification
```
for f in IPF/polyAsites/PAS/PA_counts/*.bedGraph; do mv "$f" "$(echo "$f" | sed s/_Fwd.bam_raw_PAS_Fwd.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*.bedGraph; do mv "$f" "$(echo "$f" | sed s/_Rev.bam_raw_PAS_Rev.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*FwdApa.RData; do mv "$f" "$(echo "$f" | sed s/_Fwd.bam_raw_PAS_Fwd.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*RevApa.RData; do mv "$f" "$(echo "$f" | sed s/_Rev.bam_raw_PAS_Rev.RData_/_raw_PAS_/)"; done

mv IPF/polyAsites/PAS/PA_counts/*.bedGraph IPF/polyAsites/PAS/PA_counts/bedGraphs/
```
- Sort the bedgraphs and convert them to bigWigs
```
for f in IPF/polyAsites/PAS/PA_counts/bedGraphs/*.bedGraph ; do sort -k1,1 -k2,2n "$f" -o "$f"; done
for f in IPF/polyAsites/PAS/PA_counts/bedGraphs/* ; do /home/micgdu/kentutils/bedGraphToBigWig "$f" /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome "$f.bw" ; done
mv IPF/polyAsites/PAS/PA_counts/bedGraphs/*.bw IPF/polyAsites/PAS/PA_counts/bigWigs/
for f in IPF/polyAsites/PAS/PA_counts/bigWigs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph//)"; done
```
# 9. Alternative PolyAdenylation (APA) analysis
- [7.APA.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/7.APA.R)
- APA Analysis with DEXSeq: Differential PAS usage was quantified using DEXSeq. Only active protein-coding genes (n=8629), non-overlapping on the same strand and isolated by 6kb, were analyzed
- Gene coordinates were extended 6kb downstream of the 3’ end to capture PAS beyond annotated gene ends
- Classification of APA Shifts: (DEXSeq provides log2 fold change and adjusted p-values for each PAS)
  Genes with no significant PAS change (padj ≥ 0.05) were labeled "APA no shift."
  The two most differentially used PAS were compared for genes with multiple significant PAS (padj < 0.05)
  If the distal-to-proximal PAS ratio was higher in experimental samples, the shift was distal; otherwise, it was proximal
  
# 10. Merge replicates of each cell line
- [8.Merge_Replicates.sh](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/8.Merge_Replicates.sh)







