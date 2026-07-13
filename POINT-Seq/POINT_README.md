# POINT-seq (Polymerase intact nascent transcript sequencing) Analysis
- POINT-Seq has been performed on Colorectal cells. Please refer to the script for all sample processing [1.POINT-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/POINT-Seq/1.POINT-Seq_processing.py)
- [Samples_POINT](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/POINT-Seq/Samples_POINT) - contains the name of fastq files used to process all samples at once
# Folder structure
	mkdir FastQ QC MultiQC Trimmed Aligned Replicates
	mkdir QC/FastQC QC/TrimmedQC QC/Aligned
	mkdir Aligned/ReadCounts Replicates/Merged_ReadCounts/ Replicates/Merged-BAMs 
	mkdir Replicates/Stranded_BAMs/ Replicates/Stranded_bedGraphs/ Replicates/Normalised_bedGraphs Replicates/Stranded_bigWigs/

# 1. Quality check
     fastqc Fastq/*.fq.gz -o QC/FastQC/

# 2. Adapter and low-quality reads removal
- [1.POINT-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/POINT-Seq/1.POINT-Seq_processing.py)
  
	mv Trimmed/*.zip QC/TrimmedQC/
	mv Trimmed/*.html QC/TrimmedQC/
	mv Trimmed/*.txt QC/TrimmedQC/

# 3. Alignment
- STAR aligner, human genome (38)
- Allowing for one alignment to the reference; generation of readCounts
- [1.POINT-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/POINT-Seq/1.POINT-Seq_processing.py

# 4. Merge replicates 
- [2.Merge_BAMs-Replicates.sh](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/POINT-Seq/2.Merge_BAMs-Replicates.sh)

# 5. Generate strand-specific BAMs (NEBNext Ultra II Directional RNA libraries are reverse-stranded [UTP/fr-firststrand])
- Forward transcript strand (+)
- Reads originating from transcripts on the + genomic strand are represented by: -f => includes, -F => excludes; - f 83: R1 mapped to the reverse strand and - f 163 : R2 mapped to the forward strand
- Reverse transcript strand (-)
- Reads originating from transcripts on the - genomic strand are represented by: - f 99: R1 mapped to the forward strand and - f 147 : R2 mapped to the reverse strand (-f => includes, -F => excludes)

- [3.Strand-Specific_BAMs.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/POINT-Seq/3.BAM_Reads-split.py)

# 6. Normalisation of strand specific reads
- Estimates size factors in Deseq2 and normalisation factors (per million factors)
- Generates normalised bedgraphs for full read
- [4.Normalised_bedGraphs.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/POINT-Seq/4.Normalised_bedGraphs.R)

# 7. Strand-specific bigwigs generation
- Converting normalised bedGraphs to bigwigs
	for f in Replicates/Normalised_bedGraphs/*.bedgraph; do sort -k1,1 -k2,2n "$f" -o "$f"; done

	for f in Replicates/Normalised_bedGraphs/*.bedgraph; do /home/micgdu/kentutils/bedGraphToBigWig "$f" /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome "Replicates/Stranded_bigWigs/$(basename "$f" _norm.bedgraph).bw"; done
 
# 8. Metaplots
- computeMatrix was used in combination with plotProfile as per instructions from the [deepTools](https://github.com/deeptools/deepTools)
