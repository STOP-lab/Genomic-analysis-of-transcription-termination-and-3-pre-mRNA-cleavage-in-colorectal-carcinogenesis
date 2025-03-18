# ChIP-Seq Analysis
# Folder structure
- mkdir FastQ QC MultiQC Trimmed Aligned Replicates Peaks
- mkdir QC/FastQC QC/TrimmedQC QC/Aligned
- mkdir Aligned/SplitBAM Aligned/Split_BAM/bigWigs
- mkdir Replicates/bigWigs

# 1. QUALITY CHECK
     fastqc Fastq/*.fq.gz -o QC/FastQC/
   
# 2. ADAPTER and LOW QUALITY READS REMOVAL => Script_name.py
   - mv Trimmed/*.zip QC/TrimmedQC/
   - mv Trimmed/*.html QC/TrimmedQC/
   - mv Trimmed/*.txt QC/TrimmedQC/
   
# 3. ALIGNMENT
     - #ColoRectal libraries - hybrid genome hg38_mm39 => Script_name.py
        Generated with original hg38 fasta file with chromosome naming "chr1", "chr2", "chr3" ... and modified mm39 fasta file with chromosome naming with extra "m" at the front: "mchr1", 
        "mchr2", "mchr3" ...
     - #HeLa libraries - human genome (hg38)
     - Mark duplicates
     - Index BAM file
 
# 4. SEPARATE HUMAN (experimental) AND MOUSE (spike-in) ALIGNED READS => Script_name.py
     - Mouse_Chr - contains the modified mouse chromosome names: "mchr1", "mchr2", "mchr3"... (see Mouse_Chr) which makes it easier to separate the human and mouse-aligned reads
     
# 5. GENOME COVERAGE FILES (BAM to bigWig Conversion) => Script_name.py

# 6. MERGE REPLICATES => Script_name.sh
   
# 7. PEAK CALLING WITH MACS2 => Script_name.sh

# 8. METAPLOTS
     - computeMatrix was used in combination with plotProfile as per instructions from the deepTools web page 
