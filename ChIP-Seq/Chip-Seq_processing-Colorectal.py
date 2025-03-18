#!usr/bin/env python
import os
for srr in open("Samples_Colorectal.txt"):
	sra=srr.strip("\n")
	
##Trimming adapters
	os.system("trim_galore -j 4 --paired FastQ/"+sra+"_R1.fq.gz FastQ/"+sra+"_R2.fq.gz -o Trimmed/ --basename "+sra+" --fastqc")

##Alignment 
	os.system("bowtie2 --threads 24 -x /home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2_hg38_mm39/hg38_mm39 -1 Trimmed/"+sra+"_val_1.fq.gz -2 Trimmed/"+sra+"_val_2.fq.gz --no-mixed --no-discordant 2>QC/Aligned/"+sra+"_stats.file | samtools view -Sbh | samtools sort -t 12 -o Aligned/"+sra+".bam")
 
##Mark duplicates and index the BAMs
	os.system("/home/micgdu/software/sambamba/sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/"+sra+".bam Aligned/"+sra+"_markdup.bam")
 
	os.system("samtools index -b Aligned/"+sra+"_markdup.bam") 
  
##Separate human and mouse aligned reads
#Human or experimental genome aligned reads
	os.system("samtools view -@ 15 -h Aligned/"+sra+"_markdup.bam | grep -v -f /dysk2/groupFolders/deepshika/GenomicData/Mouse/Mouse_Chr | samtools view -@ 15 -bhS - > Aligned/Split_BAM/"+sra+"-Human.bam")
	
    os.system("samtools index -@ 15 -b Aligned/Split_BAM/"+sra+"-Human.bam")

#Mouse or spike-In genome aligned reads
	os.system("samtools view -@ 15 -h Aligned/"+sra+"_markdup.bam | grep -f /dysk2/groupFolders/deepshika/GenomicData/Mouse/Mouse_Chr | samtools view -@ 15 -bhS - > Aligned/Split_BAM/"+sra+"-Mouse.bam")
	
    os.system("samtools index -@ 15 -b Aligned/Split_BAM/"+sra+"-Mouse.bam")

##BAM to bigWig Conversion
	os.system("bamCoverage -b Aligned/Split_BAM/"+sra+"-Human.bam --binSize 50 --normalizeUsing CPM -p 18 -o Aligned/bigWigs/"+sra+"-Human.bw")
