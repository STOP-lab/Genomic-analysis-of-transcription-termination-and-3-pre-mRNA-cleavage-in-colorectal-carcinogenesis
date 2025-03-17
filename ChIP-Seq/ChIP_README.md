
# Analysis folder structure
mkdir FastQ QC MultiQC Trimmed Aligned Replicates Peaks
mkdir QC/FastQC QC/TrimmedQC QC/Aligned
mkdir MultiQC/FastMQC MultiQC/TrimmedMQC MultiQC/AlignedM
mkdir Aligned/SplitBAM Aligned/Split_BAM/bigWigs
mkdir Replicates/bigWigs

 - 1. Check the quality of reads with FastQC
     fastqc Fastq/*.fq.gz -o QC/FastQC/
     multiqc  QC/FastQC/ -n RawData_multiQC -o  MultiQC/FastMQC/

2. Removal of adapters using TrimGalore

trim_galore -j 4 --paired FastQ/1CT_PCF11_rep1_R1.fq.gz FastQ/1CT_PCF11_rep1_R2.fq.gz -o Trimmed/ --basename --fastqc

mv Trimmed/*.zip QC/TrimmedQC/
mv Trimmed/*.html QC/TrimmedQC/
mv Trimmed/*.txt QC/TrimmedQC/

multiqc QC/TrimmedQC/ -n TrimmedData_multiQC -o  MultiQC/TrimmedMQC/

3. Align ChIP reads to the hybrid genome using bowtie2
### hybrid genome hg38_mm39
Made with:
- original hg38 fasta file with chromsome naming "chr1", "chr2", "chr3" ... 
- modified mm39 fasta file with chromosome naming with extra "m" at the front: "mchr1", "mchr2", "mchr3" ...

#ColoRectal libraries
bowtie2 --threads 24 -x /home/micgdu/GenomicData/genomicIndices/hsapiens/bowtie2_hg38_mm39/hg38_mm39 -1 Trimmed/1CT_PCF11_rep1_val_1.fq.gz -2 Trimmed/1CT_PCF11_rep1_val_2.fq.gz --no-mixed --no-discordant 2> QC/Aligned/1CT_PCF11_rep1_stats.txt | samtools view -Sbh | samtools sort -t 12 -o Aligned/1CT_PCF11_rep1.bam

multiqc QC/Aligned/ -n Aligned_multiQC -o  MultiQC/AlignedM/
 
#Mark duplicates using sambamba 

sambamba-1.0.1-linux-amd64-static markdup -t 25 Aligned/1CT_PCF11_rep1.bam Aligned/1CT_PCF11_rep1_markdup.bam
 
#Index BAM file using SAM tools

samtools index -@ 12 -b Aligned/1CT_PCF11_rep1_markdup.bam
  
4. Extracting the genome specific aligned reads

#Mouse_Chr - contains the modified mouse chromosome names : "mchr1", "mchr2", "mchr3" ...
Makes it easier to separate the human and mouse aligned reads

#Human or Experimental genome aligned reads
samtools view -@ 15 -h Aligned/1CT_PCF11_rep1.bam_markdup.bam | grep -v -f /dysk2/groupFolders/deepshika/GenomicData/Mouse/Mouse_Chr | samtools view -@ 15 -bhS - > Aligned/Split_BAM/1CT_PCF11_rep1-Human.bam

#Mouse or SpikeIn genome aligned reads
samtools view -@ 15 -h Aligned/1CT_PCF11_rep1_markdup.bam | grep -f /dysk2/groupFolders/deepshika/GenomicData/Mouse/Mouse_Chr | samtools view -@ 15 -bhS - > Aligned/Split_BAM/1CT_PCF11_rep1-Mouse.bam

##Indexing the BAM files
samtools index -@ 15 -b Aligned/Split_BAM/1CT_PCF11_rep1-Human.bam
  
samtools index -@ 15 -b Aligned/Split_BAM/1CT_PCF11_rep1-Mouse.bam

5. Genome Coverage files (BAM to bigWig Conversion)

bamCoverage -b Aligned/Split_BAM/1CT_PCF11_rep1-Human.bam --binSize 50 --normalizeUsing CPM -p 18 -o Aligned/Split_BAM/bigWigs/1CT_PCF11_rep1-Human.bw

6. Merge Replicates

samtools merge --threads 12 -hb Replicates/1CT_PCF11_Reps-Merged.bam Aligned/Split_BAM/1CT_PCF11_rep1-Human.bam Aligned/Split_BAM/1CT_PCF11_rep2-Human.bam Aligned/Split_BAM/1CT_PCF11_rep3-Human.bam

samtools index -@ 15 -b Replicates/1CT_PCF11_Reps-Merged.bam

bamCoverage -b Replicates/1CT_PCF11_Reps-Merged.bam --binSize 50 --normalizeUsing CPM -p 18 -o Replicates/bigWigs/1CT_PCF11_Reps-MergedBAM.bw

7. Peak calling with MACS2

macs2 callpeak -t Replicates/1CT_PCF11_Reps-Merged.bam -c Replicates/1CT_input_Reps-Merged.bam -f BAMPE -g hs --outdir  Peaks/ -n 1CT_PCF11_Reps-MergedBAM --nomodel --broad -q 0.01 --broad-cutoff 0.01 2> Peaks/1CT_PCF11_Reps-MergedBAM_macs2.log

8. For Metaplots and heatmaps, computeMatrix was used in combination with plotProfile/plotHeatmap as per instructions from deepTools web page


