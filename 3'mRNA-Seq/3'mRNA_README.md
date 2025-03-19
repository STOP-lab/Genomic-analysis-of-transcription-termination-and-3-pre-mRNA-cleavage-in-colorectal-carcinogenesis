# 3'mRNA-Seq Analysis
- mNET-Seq has been performed on Colorectal cells. Please refer to the script for all sample processing [1.3Prime-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/1.3Prime-Seq_processing.py)
- [Samples_3'mRNA](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/Samples_3'mRNA) - contains the name of fastq files used to process all samples at once
# Folder structure for analysis
	mkdir 3Prime-Seq/FastQ 3Prime-Seq/QC 3Prime-Seq/Trimmed 3Prime-Seq/Trimmed/PolyA-T/ 3Prime-Seq/Aligned mkdir 3Prime-Seq/Aligned/ReadCounts/
 	mkdir 3Prime-Seq/QC/FastQC 3Prime-Seq/QC/TrimmedQC 3Prime-Seq/QC/TrimmedQC/PolyA-T 3Prime-Seq/QC/Aligned
	mkdir 3Prime-Seq/Replicates 3Prime-Seq/Aligned/Stranded_BAMs/ 3Prime-Seq/Aligned/Stranded_bedGraphs/ 
	mkdir 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs 3Prime-Seq/Aligned/Stranded_bigWigs  
# 1. DeMultiplexing and FastQ format conversion
  	bcl2fastq -R 231030_A01680_0155_AHMN5MDRX3/ -o 231030_A01680_0155_AHMN5MDRX3/3Prime-Seq/FastQ/ --no-lane-splitting -p 26
- Raw .cbcl files to fastq conversion and demultiplexing  was done using custom-created SampleSheet
- --no-lane-splitting Make sure that the reads of one sample from different lanes don't produce multiple fastq files but would end up in one FastQ file
# 2. Quality check
	fastqc 3Prime-Seq/FastQ/*.gz -o 3Prime-Seq/QC/FastQC/
# 3. Adapter, low-quality reads, and PolyA/T_tails removal
	mv Trimmed/*.zip QC/TrimmedQC/
	mv Trimmed/*.html QC/TrimmedQC/
	mv Trimmed/*.txt QC/TrimmedQC/
	mv Trimmed/*.zip QC/TrimmedQC/PolyA-T
	mv Trimmed/*.html QC/TrimmedQC/PolyA-T
	mv Trimmed/*.txt QC/TrimmedQC/PolyA-T
 - Move the trimmed reports to the quality check reports to main the folder structure
 - [1.3Prime-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/3'mRNA-Seq/1.3Prime-Seq_processing.py)
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
- Sort normalised bedgraphs and generate bigWigs
```
for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do sort -k1,1 -k2,2n "$f" -o "$f"; done
for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do /home/micgdu/kentutils/bedGraphToBigWig "$f" hg38_chromsizes.genome "$f.bw"; done
for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph.bw/.bw/g)"; done
mv 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw 3Prime-Seq/Aligned/Stranded_bigWigs/

6. Filtering Internal Primming Events => InternalPrimming_Mask.sh; InternalPrimming_Filtering.R

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do sort -k1,1 -k2,2n "$f" -o "$f"; done

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do /home/micgdu/kentutils/bedGraphToBigWig "$f" hg38_chromsizes.genome "$f.bw"; done

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph.bw/.bw/g)"; done

mv 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw 3Prime-Seq/Aligned/Stranded_bigWigs/

7. Polyadenylation site (PAS) Detection and Quantification => PAS_Calling.R; PAS_Quantification.R

######################## Renaming and relocating files after PAS Quantification ########################
for f in IPF/polyAsites/PAS/PA_counts/*.bedGraph; do mv "$f" "$(echo "$f" | sed s/_Fwd.bam_raw_PAS_Fwd.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*.bedGraph; do mv "$f" "$(echo "$f" | sed s/_Rev.bam_raw_PAS_Rev.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*FwdApa.RData; do mv "$f" "$(echo "$f" | sed s/_Fwd.bam_raw_PAS_Fwd.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*RevApa.RData; do mv "$f" "$(echo "$f" | sed s/_Rev.bam_raw_PAS_Rev.RData_/_raw_PAS_/)"; done

mv IPF/polyAsites/PAS/PA_counts/*.bedGraph IPF/polyAsites/PAS/PA_counts/bedGraphs/

######### Sorting the bedgraphs and converting them to bigWigs
for f in IPF/polyAsites/PAS/PA_counts/bedGraphs/*.bedGraph ; do sort -k1,1 -k2,2n "$f" -o "$f"; done
for f in IPF/polyAsites/PAS/PA_counts/bedGraphs/* ; do /home/micgdu/kentutils/bedGraphToBigWig "$f" /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome "$f.bw" ; done
mv IPF/polyAsites/PAS/PA_counts/bedGraphs/*.bw IPF/polyAsites/PAS/PA_counts/bigWigs/
for f in IPF/polyAsites/PAS/PA_counts/bigWigs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph//)"; done

8. Alternative PolyAdenylation (APA) analysis => APA.R

9. Merge replicates of each cell line => Merge_BAMs-Replicates.sh

10. For Metaplots and heatmaps, computeMatrix was used in combination with plotProfile/plotHeatmap as per instructions from deepTools web page





