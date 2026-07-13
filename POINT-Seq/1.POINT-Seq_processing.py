#!usr/bin/env python
import os
for srr in open("Samples_POINT.txt"):
	sra=srr.strip("\n")

### Adapter trimming : TrimGalore
	os.system("trim_galore -j 4 --paired FastQ/"+sra+"_R1.fq.gz FastQ/"+sra+"_R2.fq.gz --basename "+sra+" --fastqc_args \"-t 18 --outdir QC/Trimmed_FastQC/\" -o Trimmed_FastQ/")

### Alignment: STAR
	os.system("STAR --genomeDir /home/micgdu/GenomicData/genomicIndices/hsapiens/STAR/STAR-2.7.6a_GRCh38_120bp/ --runThreadN 10 --readFilesCommand zcat --readFilesIn Trimmed_FastQ/"+sra+"_val_1.fq.gz Trimmed_FastQ/"+sra+"_val_2.fq.gz --outFileNamePrefix Aligned/"+sra+"-STAR --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD --outMultimapperOrder Random --outSAMmultNmax 1 --chimSegmentMin 15 --quantMode GeneCounts")
 
### Mark duplicates: sambamba
	os.system("/home/micgdu/software/sambamba/sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/"+sra+"-STARAligned.sortedByCoord.out.bam Aligned/"+sra+"_markdup.bam")

### Index BAM: samtools
	os.system("samtools index -@ 12 -b Aligned/"+sra+"_markdup.bam")