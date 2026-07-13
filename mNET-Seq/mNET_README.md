# mNET-Seq (Mammalian Native elongating transcript sequencing) Analysis
- mNET-Seq has been performed on Colorectal and Pancreatic cells. For Pancreatic cells, only read2 was considered for downstream analysis because read1 was of poor technical quality.
- Please refer to the scripts in [mNET-Seq](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/) for complete analysis
- [Samples_mNET](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/Samples_mNET) - contains the names of the fastq files used to process all samples at once
# Folder structure
	mkdir FastQ QC MultiQC Trimmed Aligned Replicates
	mkdir QC/FastQC QC/TrimmedQC QC/Aligned
	mkdir Aligned/ReadCounts Replicates/Merged_ReadCounts/ Replicates/Merged-BAMs 
	mkdir Replicates/Stranded_BAMs/ Replicates/Stranded_bedGraphs/ Replicates/Stranded_bedGraphs/Normalised_bedGraphs
	mkdir Replicates/SNR/Stranded_bedGraphs/ Replicates/SNR/Normalised_bedGraphs mkdir Replicates/SNR/Stranded_bigWigs/
	
# 1. Quality check
     fastqc Fastq/*.fq.gz -o QC/FastQC/
# Raw FASTQ stats
	multiqc QC/FastQC/ -o MultiQC -n Raw_multiqc
	
# 2. Adapter and low-quality reads removal
- [1.mNET-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/1.mNET-Seq_processing.py)
 # Trimmed FASTQ stats 
	mv Trimmed/*.zip QC/TrimmedQC/
	mv Trimmed/*.html QC/TrimmedQC/
	mv Trimmed/*.txt QC/TrimmedQC/
	multiqc QC/TrimmedQC/ -o MultiQC -n Trimmed_multiqc
	
# 3. Alignment
- STAR aligner, human genome (38)
- Allowing for one alignment to the reference; generation of read counts
- [1.mNET-Seq_processing.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/1.mNET-Seq_processing.py)
# Alignment Stats
	multiqc Aligned/*Log.final.out -o MultiQC -n Aligned_multiqc
	
# 4. Merge replicates 
- [2.Merge_BAMs-Replicates.sh](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/2.Merge_BAMs-Replicates.sh)
# For Pancreatic cells, use the commands below
	samtools merge -@ 18 -c Replicates/Merged_BAMs/MiaPaCa2-T4ph-SE.bam Aligned/MiaPaCa2-T4ph_rep3-SE.bam Aligned/MiaPaCa2-T4ph_rep4-SE.bam
	samtools merge -@ 18 -c Replicates/Merged_BAMs/Panc1-T4ph-SE.bam Aligned/Panc1-T4ph_rep3-SE.bam Aligned/Panc1-T4ph_rep4-SE.bam
	for f in Replicates/Merged_BAMs/*.bam; do samtools index -@ 12 -b "$f"; done;

# 5. Generate strand-specific BAMs
- Forward Strand reads - 5' to 3' direction
	-f => includes, -F => excludes; 64 => First in pair (Forward Strand reads - 5' to 3' direction), 128 => Second in pair, 16 => read reverse strand, 144= 128+16
- Reverse Strand reads - 3' to 5' direction)
	-f => includes, -F => excludes; 128 => Second in pair (Reverse Strand reads - 3' to 5' direction), 16 => read reverse strand, 80 => First in pair and read reverse strand
- [3.BAM_Reads-split.py](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/3.BAM_Reads-split.py)

# For Pancreatic cells, use the commands below
	for f in Replicates/Merged_BAMs/*.bam; do base=$(basename "$f" .bam); samtools view -bh -@ 18 -f 16 -F 4 "$f" > Replicates/Stranded_BAMs/"${base}_fwd.bam"; done
	for f in Replicates/Merged_BAMs/*.bam; do base=$(basename "$f" .bam); samtools view -bh -@ 18 -F 20 "$f" > Replicates/Stranded_BAMs/"${base}_rev.bam"; done
for f in Replicates/Merged_BAMs/Stranded_BAMs/*.bam; do samtools index -@ 18 -b "$f"; done;

## Generate coverage files for pancreatic samples
	for f in Replicates//Stranded_BAMs/*_fwd.bam; do f_base=$(basename "$f" .bam); bedtools genomecov -ibam  Replicates/Stranded_BAMs/"${f_base}.bam" -bga -split | gzip - > Replicates/Stranded_bedGraphs/"${f_base}.bedgraph.gz"; done
	for f in Replicates/Stranded_BAMs/*_rev.bam; do f_base=$(basename "$f" .bam); bedtools genomecov -ibam  Replicates/Stranded_BAMs/"${f_base}.bam" -bga -split | gzip - > Replicates/Stranded_bedGraphs/"${f_base}.bedgraph.gz"; done

# 6. Normalisation of strand-specific reads
- Estimates size factors in Deseq2 and normalisation factors (per million factors)
- Generates normalised bedgraphs for full read
- [4.Normalised_bedGraphs.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/4.Normalised_bedGraphs.R)
- For Pancreatic cells [4.Normalised_bedGraphs_Pancreatic.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/4.Normalised_bedGraphs_Pancreatic.R)

# Merge the readcounts of replicates before normalisation
	paste Aligned/ReadCounts/MiaPaCa2-T4ph_rep1-STAR-SEReadsPerGene.out.tab Aligned/ReadCounts/MiaPaCa2-T4ph_rep2-STAR-SEReadsPerGene.out.tab | awk '{printf "%s", $1; for (c=2; c<=4; c++) {s=0; for(i=c;i<=NF;i+=4)s+=$i; printf "\t%s", s;} print ""}' > Replicates/Merged_ReadCounts/MiaPaCa2-T4ph-STAR-SEReadsPerGene.out.tab
	paste Aligned/ReadCounts/Panc1-T4ph_rep1-STAR-SEReadsPerGene.out.tab Aligned/ReadCounts/Panc1-T4ph_rep2-STAR-SEReadsPerGene.out.tab | awk '{printf "%s", $1; for (c=2; c<=4; c++) {s=0; for(i=c;i<=NF;i+=4)s+=$i; printf "\t%s", s;} print ""}' > Replicates/Merged_ReadCounts/Panc1-T4ph-STAR-SEReadsPerGene.out.tab

# 7. Extract the single nucleotide resolution (SNR)
- Extracts the last transcribed nucleotide from a read, that is, the last position of the read
- [5.SNR.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/5.SNR.R)

# 8. Normalisation of SNR reads
- Estimates size factors in Deseq2 and normalisation factors (per million factors)
- Generates normalised bedgraphs for single nucleotide resolution reads
- [6.SNR-Normalised_bedGraphs.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/6.SNR-Normalised_bedGraphs.R)
- For Pancreatic cells [6.SNR-Normalised_bedGraphs_Pancreatic.R](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/6.SNR-Normalised_bedGraphs_Pancreatic.R)
  
# 9. Termination Windows
- MACS2 bdgbroadcall function was called upon T4ph stranded bedgraphs with default paprameters
- Reverse strand bedgraph files were multiplied by -1 as they contain negative scores, and that might interfere with identifying the significantly enriched T4ph regions
- [7.Termination_Windows.sh](https://github.com/STOP-lab/Genomic-analysis-of-transcription-termination-and-3-pre-mRNA-cleavage-in-colorectal-carcinogenesis/blob/main/mNET-Seq/7.Termination_Windows.sh)
  
# 10. Metaplots
- computeMatrix was used in combination with plotProfile as per instructions from the [deepTools](https://github.com/deeptools/deepTools)
