#!usr/bin/env python
import os
for srr in open("Samples_HeLa.txt"):
	sra=srr.strip("\n")
	
##Trimming adapters
	os.system("trim_galore -j 4 --paired FastQ/"+sra+"_R1.fastq.gz FastQ/"+sra+"_R2.fastq.gz -o Trimmed/ --basename "+sra+" --fastqc")

##Alignment 
 	os.system("bowtie2 --threads 24 -x /home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2/hg38 -1 Trimmed/"+sra+"_val_1.fq.gz -2 Trimmed/"+sra+"_val_2.fq.gz --no-mixed --no-discordant 2>QC/Aligned/"+sra+"_stats.file | samtools view -Sbh | samtools sort -t 12 -o Aligned/"+sra+".bam")
 
##Mark duplicates and index the BAMs
	os.system("/home/micgdu/software/sambamba/sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/"+sra+".bam Aligned/"+sra+"_markdup.bam")
 
	os.system("samtools index -b Aligned/"+sra+"_markdup.bam") 
  
##BAM to bigWig Conversion
	os.system("bamCoverage -b Aligned/Split_BAM/"+sra+"_markdup.bam --binSize 50 --normalizeUsing CPM -p 18 -o Aligned/bigWigs/"+sra+".bw")