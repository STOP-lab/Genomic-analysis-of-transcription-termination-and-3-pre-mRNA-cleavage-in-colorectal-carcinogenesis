#!usr/bin/env python
import os
for srr in open("Samples_3mRNA.txt"):
	sra=srr.strip("\n")

### Adapter and PolyA/T_tail trimming using cutadapt
	os.system("cutadapt -j 4 -q 20 -m 20 -O 10 --times 2 --nextseq-trim=10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a A{20} -g T{20} -o Trimmed/"+sra+"_trimmed_clean.fq.gz FastQ/"+sra+"_R1_001.fastq.gz > Trimmed/"+sra+"_cutadapt_report.txt")

## Alignment using STAR 
	os.system("STAR --genomeDir /home/micgdu/GenomicData/genomicIndices/hsapiens/STAR/STAR-2.7.6a_GRCh38_120bp/ --runThreadN 15 --readFilesCommand zcat --readFilesIn Trimmed/PolyA/"+sra+"_trimmed_clean.fq.gz --outFileNamePrefix Aligned/"+sra+"-STAR --outFilterType BySJout --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD --quantMode GeneCounts")

## Mark the duplicates using sambamba
	os.system("/home/micgdu/software/sambamba/sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/"+sra+"-STARAligned.sortedByCoord.out.bam Aligned/"+sra+"_markdup.bam")
	os.system("samtools index -@ 12 -b Aligned/"+sra+"_markdup.bam")

## Generate Strand-specfic BAMs
	os.system("samtools view -bh -@ 12 -f 16 Aligned/"+sra+"_markdup.bam -o Aligned/Stranded_BAMs/"+sra+"_Fwd.bam")
	os.system("samtools view -bh -@ 12 -F 16 Aligned/"+sra+"_markdup.bam -o Aligned/Stranded_BAMs/"+sra+"_Rev.bam")
	os.system("samtools index -@ 12 -b Aligned/Stranded_BAMs/"+sra+"_Fwd.bam")
	os.system("samtools index -@ 12 -b Aligned/Stranded_BAMs/"+sra+"_Rev.bam")

## Genome Coverage files (BAM to bedGraph Conversion)
	os.system("bedtools genomecov -ibam Aligned/Stranded_BAMs/"+sra+"_Fwd.bam -bga -split | gzip - > Aligned/Stranded_bedGraphs/"+sra+"_Fwd.bedgraph.gz")
	os.system("bedtools genomecov -ibam Aligned/Stranded_BAMs/"+sra+"_Rev.bam -bga -split | gzip - > Aligned/Stranded_bedGraphs/"+sra+"_Rev.bedgraph.gz")
