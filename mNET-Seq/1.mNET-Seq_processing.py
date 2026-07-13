#!usr/bin/env python
import os
for srr in open("Samples_mNET.txt"):
	sra=srr.strip("\n")

################### Colorectal mNET-Seq Samples	(Paired-End mode) ###################
### Trimming of adapters
	os.system("trim_galore -j 4 --paired FastQ/"+sra+"_R1.fastq.gz FastQ/"+sra+"_R2.fastq.gz -o Trimmed/ --basename "+sra+" --fastqc")

### Alignment
	os.system("STAR --genomeDir /home/micgdu/GenomicData/genomicIndices/hsapiens/STAR/STAR-2.7.6a_GRCh38_100bp/ --runThreadN 10 --readFilesCommand zcat --readFilesIn Trimmed/"+sra+"_val_1.fq.gz Trimmed/"+sra+"_val_2.fq.gz --outFileNamePrefix Aligned/"+sra+"-STAR --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD --outMultimapperOrder Random --outSAMmultNmax 1 --chimSegmentMin 15 --quantMode GeneCounts")

### Mark duplicates using sambamba
	os.system("/home/micgdu/software/sambamba/sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/"+sra+"-STARAligned.sortedByCoord.out.bam Aligned/Markdup/"+sra+"_markdup.bam")	
	 
##Index BAM file using SAM tools
	os.system("samtools index -@ 12 -b Aligned/"+sra+"-STARAligned.sortedByCoord.out.bam")
	
################### Pancreatic T4ph mNET-Seq Samples (Single-End (SE) Mode was chosen as R1 quality was technically poor) ###################
	os.system("trim_galore -j 4 --nextseq 10 FastQ/"+sra+"_R2.fq.gz -o  Trimmed/ --basename "+sra+" --fastqc_args \"--outdir QC/Trimmed_FastQC/\" --adapter GATCGTCGGACTGTAGAACTCTGAAC")
 
### Alignment
	os.system("STAR --runThreadN 12 --genomeDir /home/micgdu/GenomicData/genomicIndices/hsapiens/STAR/STAR-2.7.6a_GRCh38_120bp/ --readFilesIn Trimmed/"+sra+"_trimmed.fq.gz  --readFilesCommand zcat --outFileNamePrefix Aligned/"+sra+"-STAR-SE --outSAMtype BAM SortedByCoordinate --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMismatchNmax 2 --outFilterMatchNmin 15 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 -outSAMattributes NH HI NM MD --outMultimapperOrder Random --outSAMmultNmax 1 --chimSegmentMin 15 --quantMode GeneCounts")

### Mark duplicates using sambamba
	os.system("/home/micgdu/software/sambamba/sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/"+sra+"-STAR-SEAligned.sortedByCoord.out.bam Aligned/Markdup/"+sra+"-SE_markdup.bam")

##Index BAM file using SAM tools
	os.system("samtools index -@ 12 -b Aligned/Markdup/"+sra+"-SE_markdup.bam"")

 
