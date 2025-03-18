# mNET-Seq Analysis
- mNET-Seq has been performed on Colorectal cells. Please refer to the script for all sample processing [Chip-Seq_processing-Colorectal.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/ChIP-Seq/Chip-Seq_processing-Colorectal.py)

- [Samples_Colorectal.txt](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/ChIP-Seq/Samples_Colorectal.txt) - contains the name of fastq files used to process all samples at once
# Folder structure
	mkdir FastQ QC MultiQC Trimmed Aligned Replicates
	mkdir QC/FastQC QC/TrimmedQC QC/Aligned
	mkdir Aligned/ReadCounts Replicates/Merged_ReadCounts/ Replicates/Merged-BAMs 
	mkdir Replicates/Stranded_BAMs/ Replicates/Stranded_bedGraphs/ Replicates/Stranded_bedGraphs/Normalised_bedGraphs
	mkdir Replicates/Stranded_bedGraphs/SNR/ Replicates/Stranded_bedGraphs/SNR/Normalised_bedGraphs
# 1. Quality check
     fastqc Fastq/*.fq.gz -o QC/FastQC/
# 2. Adapter and low-quality reads removal
	mv Trimmed/*.zip QC/TrimmedQC/
	mv Trimmed/*.html QC/TrimmedQC/
	mv Trimmed/*.txt QC/TrimmedQC/
- [Chip-Seq_processing-Colorectal.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/ChIP-Seq/Chip-Seq_processing-Colorectal.py)
# 3. Alignment
- STAR aligner, human genome (38)
- Allowing for one alignment to the reference; generation of readCounts[Chip-Seq_processing-HeLa.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/ChIP-Seq/Chip-Seq_processing-HeLa.py)
# 4. Merge replicates 
- [MergeRepsnPeakCalling.sh](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/ChIP-Seq/MergeRepsnPeakCalling.sh)
# 5. Generate strand-specific BAMs
- Forward Strand reads - 5' to 3' direction
	-f => includes, -F => excludes; 64 => First in pair (Forward Strand reads - 5' to 3' direction), 128 => Second in pair, 16 => read reverse strand, 144= 128+16
- Reverse Strand reads - 3' to 5' direction)
	-f => includes, -F => excludes; 128 => Second in pair (Reverse Strand reads - 3' to 5' direction), 16 => read reverse strand, 80 => First in pair and read reverse strand
- [Chip-Seq_processing-Colorectal.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/ChIP-Seq/Chip-Seq_processing-Colorectal.py)
# 6. Normalisation of strand specific reads
- Estimates size factors in Deseq2 and normalisation factors (per million factors)
- Generates normalised bedgraphs for full read
# 7. Extracting the single nucelotide resolution (SNR)
- Extracts the last transcribed nucleotide from a read, that is the last position of read
- [MergeRepsnPeakCalling.sh](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/ChIP-Seq/MergeRepsnPeakCalling.sh)
# 8. Normalisation of SNR reads
- Estimates size factors in Deseq2 and normalisation factors (per million factors)
- Generates normalised bedgraphs for single nucleotide resoultion reads
# 9. Termination Windows
- MACS2 bdgbroadcall function was called upon T4ph stranded bedgraphs with default paprameters
- Reverse strand bedgraph files were multiplied with -1 as the contain negative scores and that might interfere in identifying the significantly enriched T4ph regions
# 10. Metaplots
- computeMatrix was used in combination with plotProfile as per instructions from the [deepTools](https://github.com/deeptools/deepTools)
