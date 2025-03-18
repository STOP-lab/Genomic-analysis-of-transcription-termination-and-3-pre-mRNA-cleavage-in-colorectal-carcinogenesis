# ChIP-Seq Analysis
# Folder structure
     mkdir FastQ QC MultiQC Trimmed Aligned Replicates Peaks
     mkdir QC/FastQC QC/TrimmedQC QC/Aligned
     mkdir Aligned/SplitBAM Aligned/Split_BAM/bigWigs
     mkdir Replicates/Merged-BAMs Replicates/Merged-BAMs/bigWigs
# 1. QUALITY CHECK
     fastqc Fastq/*.fq.gz -o QC/FastQC/
# 2. ADAPTER and LOW QUALITY READS REMOVAL
     trim_galore -j 4 --paired FastQ/1CT_PCF11_rep1_R1.fq.gz FastQ/1CT_PCF11_rep1_R2.fq.gz -o Trimmed/ --basename 1CT_PCF11 --fastqc
     mv Trimmed/*.zip QC/TrimmedQC/
     mv Trimmed/*.html QC/TrimmedQC/
     mv Trimmed/*.txt QC/TrimmedQC/
# 3. ALIGNMENT
- ColoRectal libraries - hybrid genome hg38_mm39; generated with original hg38 fasta file with chromosome naming "chr1", "chr2", "chr3" ... and modified mm39 fasta file with chromosome naming with extra "m" at the front: "mchr1", "mchr2", "mchr3"...
- HeLa libraries - human genome (hg38)
# 4. SEPARATE HUMAN (experimental) AND MOUSE (spike-in) ALIGNED READS => Script_name.py
- Mouse_Chr - contains the modified mouse chromosome names: "mchr1", "mchr2", "mchr3"... (see Mouse_Chr) which makes it easier to separate the human and mouse-aligned reads
# 5. GENOME COVERAGE FILES (BAM to bigWig Conversion)
     bamCoverage -b Aligned/Split_BAM/1CT_PCF11_rep1-Human.bam --binSize 50 --normalizeUsing CPM -p 18 -o Aligned/Split_BAM/bigWigs/1CT_PCF11_rep1-Human.bw
     bamCoverage -b Aligned/HeLa_PCF11_siPCF11_markdup.bam --binSize 50 --normalizeUsing CPM -p 18 -o Aligned/bigWigs/HeLa_PCF11_siPCF11.bw
# 6. MERGE REPLICATES => Script_name.sh
     samtools merge --threads 12 -hb Replicates/Merged-BAMs/1CT_PCF11_Reps-Merged.bam Aligned/Split_BAM/1CT_PCF11_rep1-Human.bam Aligned/Split_BAM/1CT_PCF11_rep2-Human.bam Aligned/Split_BAM/1CT_PCF11_rep3-Human.bam
     samtools index -@ 15 -b Replicates/1CT_PCF11_Reps-Merged.bam
# 7. PEAK CALLING WITH MACS2 => Script_name.sh
     macs2 callpeak -t Replicates/Merged-BAMs/1CT_PCF11_Reps-Merged.bam -c Replicates/Merged-BAMs/1CT_input_Reps-Merged.bam -f BAMPE -g hs --outdir Peaks/ -n 1CT_PCF11_Reps-MergedBAM --nomodel --broad -q 0.01 --broad-cutoff 0.01 2> Peaks/1CT_PCF11_Reps-MergedBAM_macs2.log
# 8. METAPLOTS
- computeMatrix was used in combination with plotProfile as per instructions from the deepTools web page 
