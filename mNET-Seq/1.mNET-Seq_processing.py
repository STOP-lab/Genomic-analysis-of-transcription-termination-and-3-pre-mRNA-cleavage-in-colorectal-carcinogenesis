#!usr/bin/env python
import os
for srr in open("Samples_mNET.txt"):
	sra=srr.strip("\n")

### Trimming of adapters
	os.system("trim_galore -j 4 --paired FastQ/"+sra+"_R1.fastq.gz FastQ/"+sra+"_R2.fastq.gz -o Trimmed/ --basename "+sra+" --fastqc")

### Alignment
	os.system("STAR --genomeDir /home/micgdu/GenomicData/genomicIndices/hsapiens/STAR/STAR-2.7.6a_GRCh38_100bp/ --runThreadN 10 --readFilesCommand zcat --readFilesIn Trimmed/"+sra+"_val_1.fq.gz Trimmed/mNET/"+sra+"_val_2.fq.gz --outFileNamePrefix Aligned/100bp/"+sra+"-STAR --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD --outMultimapperOrder Random --outSAMmultNmax 1 --chimSegmentMin 15 --quantMode GeneCounts")
 
##Index BAM file using SAM tools
	os.system("samtools index -@ 12 -b Aligned/"+sra+"-STARAligned.sortedByCoord.out.bam")
