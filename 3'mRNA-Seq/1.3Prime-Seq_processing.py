#!usr/bin/env python
import os
for srr in open("accession_id.txt"):
	sra=srr.strip("\n")

### Adapter and PolyA/T_tail trimming using TrimGalore

##Trimming adapters
	os.system("trim_galore -j 4  3Prime-Seq/FastQ/"+sra+".fastq.gz  -o 3Prime-Seq/Trimmed/ -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --basename "+sra+" --fastqc_args '--threads 12 --outdir 3Prime-Seq/QC/Trimmed_FastQC/ ")

## Trimming Poly-A tail
	os.system("cutadapt -j 4 -m 20 -o 20 -g "A{20}" -o 3Prime-Seq/Trimmed/PolyA/"+sra+"-PolyA_CT.fq.gz 3Prime-Seq/Trimmed/"+sra+"_trimmed.fq.gz > 3Prime-Seq/Trimmed/PolyA/"+sra+"-PolyA_CT_trimming_report.txt")

## Alignment using STAR 
	os.system("STAR --genomeDir /home/micgdu/GenomicData/genomicIndices/hsapiens/STAR/STAR-2.7.6a_GRCh38_120bp/ --runThreadN 15 --readFilesCommand zcat --readFilesIn 3Prime-Seq/Trimmed/PolyA/"+sra+"-PolyA_CT.fq.gz --outFileNamePrefix 3Prime-Seq/Aligned/"+sra+"-STAR --outFilterType BySJout --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD --quantMode GeneCounts --genomeLoad LoadAndKeep")

## Mark the duplicates using sambamba
	os.system("/home/micgdu/software/sambamba/sambamba-1.0.1-linux-amd64-static markdup -t 25 3Prime-Seq/Aligned/"+sra+"-STARAligned.sortedByCoord.out.bam 3Prime-Seq/Aligned/"+sra+"_markdup.bam")
	os.system("samtools index -@ 12 -b 3Prime-Seq/Aligned/"+sra+"_markdup.bam")

## Generate Strand-specfic BAMs
	os.system("samtools view -bh -@ 12 -f 16 3Prime-Seq/Aligned/"+sra+"_markdup.bam -o 3Prime-Seq/Aligned/Stranded_BAMs/"+sra+"_Fwd.bam")
	os.system("samtools view -bh -@ 12 -F 16 3Prime-Seq/Aligned/"+sra+"_markdup.bam -o 3Prime-Seq/Aligned/Stranded_BAMs/"+sra+"_Rev.bam")
	os.system("samtools index -@ 12 -b ColoRectal_Cancer/3Prime-Seq/Aligned/Stranded_BAMs/"+sra+"_Fwd.bam")
	os.system("samtools index -@ 12 -b ColoRectal_Cancer/3Prime-Seq/Aligned/Stranded_BAMs/"+sra+"_Rev.bam")

## Genome Coverage files (BAM to bedGraph Conversion)
	os.system("bedtools genomecov -ibam 3Prime-Seq/Aligned/Stranded_BAMs/"+sra+"_Fwd.bam -bga -split | gzip - > 3Prime-Seq/Aligned/Stranded_bedGraphs/"+sra+"_Fwd.bedgraph.gz")
	os.system("bedtools genomecov -ibam 3Prime-Seq/Aligned/Stranded_BAMs/"+sra+"_Rev.bam -bga -split | gzip - > 3Prime-Seq/Aligned/Stranded_bedGraphs/"+sra+"_Rev.bedgraph.gz")
