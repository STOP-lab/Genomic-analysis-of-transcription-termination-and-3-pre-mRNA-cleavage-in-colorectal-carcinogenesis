##Analysis folder structure

mkdir FastQ QC MultiQC Trimmed Aligned Replicates Peaks
mkdir QC/FastQC QC/TrimmedQC/ QC/Aligned
mkdir MultiQC/FastMQC MultiQC/TrimmedMQC/ MultiQC/AlignedM

1. Check the quality of reads with FastQC

fastqc Fastq/*.fastq.gz -o QC/FastQC/

multiqc  QC/FastQC/ -n RawData_multiQC -o  MultiQC/FastMQC/

2. Trimming adapters

trim_galore -j 4 --paired FastQ/HeLa_PCF11_siPCF11_R1.fastq.gz FastQ/HeLa_PCF11_siPCF11_R2.fastq.gz -o Trimmed/ --basename --fastqc

mv Trimmed/*.zip QC/TrimmedQC/
mv Trimmed/*.html QC/TrimmedQC/
mv Trimmed/*.txt QC/TrimmedQC/

multiqc QC/TrimmedQC/ -n TrimmedData_multiQC -o  MultiQC/TrimmedMQC/

3. Alignment

#HeLa libraries
bowtie2 --threads 24 -x /home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2/hg38 -1 Trimmed/HeLa_PCF11_siPCF11_val_1.fq.gz -2 Trimmed/HeLa_PCF11_siPCF11_val_2.fq.gz --no-mixed --no-discordant 2>QC/Aligned/HeLa_PCF11_siPCF11_stats.txt | samtools view -Sbh | samtools sort -t 12 -o Aligned/HeLa_PCF11_siPCF11.bam

multiqc QC/Aligned/ -n Aligned_multiQC -o  MultiQC/AlignedM/
 
#Mark duplicates using sambamba 

sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/HeLa_PCF11_siPCF11.bam Aligned/HeLa_PCF11_siPCF11_markdup.bam
 
#Index BAM file using SAM tools

samtools index -@ 12 -b Aligned/HeLa_PCF11_siPCF11_markdup.bam

##Indexing the BAM files
samtools index -@ 15 -b Aligned/HeLa_PCF11_siPCF11_markdup.bam
  
5. Genome Coverage files (BAM to bigWig Conversion)

bamCoverage -b Aligned/HeLa_PCF11_siPCF11_markdup.bam --binSize 50 --normalizeUsing CPM -p 18 -o Aligned/bigWigs/HeLa_PCF11_siPCF11.bw

6. For Metaplots and heatmaps, computeMatrix was used in combination with plotProfile/plotHeatmap as per instructions from deepTools web page


