##Analysis folder structure

mkdir FastQ QC MultiQC Trimmed Aligned Replicates 
mkdir QC/FastQC QC/TrimmedQC QC/Aligned
mkdir MultiQC/FastMQC MultiQC/TrimmedMQC MultiQC/AlignedM
mkdir Aligned/ReadCounts/
mkdir Replicates/Merged_ReadCounts/ Replicates/Stranded_BAMs/ Replicates/Stranded_bedGraphs/ Replicates/Stranded_bedGraphs/Normalised_bedGraphs Replicates/Stranded_bedGraphs/SNR/ Replicates/Stranded_bedGraphs/SNR/Normalised_bedGraphs 

1. Check the quality of reads with FastQC

fastqc FastQ/*.gz -o QC/FastQC/

multiqc  QC/FastQC/ -n RawData_multiQC -o  MultiQC/FastMQC/

2. Removal of adapters using TrimGalore => 1.mNET-Seq_procesing.py

mv Trimmed/*.zip QC/TrimmedQC/
mv Trimmed/*.html QC/TrimmedQC/
mv Trimmed/*.txt QC/TrimmedQC/

multiqc QC/TrimmedQC/ -n TrimmedData_multiQC -o  MultiQC/TrimmedMQC/

3. Aligning the mNET-Seq reads to the hg38 reference genome using STAR => 1.mNET-Seq_procesing.py

mv Aligned/STARReadsPerGene.out.tab Aligned/ReadCounts/
mv Aligned/STARLog.final.out QC/Aligned

multiqc QC/Aligned/ -n Aligned_multiQC -o  MultiQC/AlignedM/

4. Merge replicates of each cell line => 2.Merge_BAMs-Replicates.sh

5. Create strand-specfic BAMs and genome coverage files(bedGraps) => 3.BAM_Reads-split.py

### Forward Strand reads - 5' to 3' direction
### -f => includes, -F => excludes; 64 => First in pair (Forward Strand reads - 5' to 3' direction), 128 => Second in pair, 16 => read reverse strand, 144= 128+16

### Reverse Strand reads - 3' to 5' direction)
### -f => includes, -F => excludes; 128 => Second in pair (Reverse Strand reads - 3' to 5' direction), 16 => read reverse strand, 80 => First in pair and read reverse strand

rm Aligned/Stranded_BAMs/*1.bam
rm Aligned/Stranded_BAMs/*2.bam

6. Full read Normalisation => 4.Normalised_bedGraphs.R

7. Single nucleotide resolution (SNR); Extracting the last transcribed nucleotide => 5.SNR.R

8. SNR_Normalisation => 6.SNR-Normalised_bedGraphs.R

mv Replicates/Stranded_bedGraphs/SNR/*_norm.bedgraph Replicates/Stranded_bedGraphs/SNR/Normalised_bedGraphs/

9. Termination Windows => 7. Termination_Windows.sh

## MACS2 bdgbroadcall fucntion was called upon stranded bedgraph files
## Reverse strand bedgraph files were multiplied with -1 as the have negative scores and that can interfere in finding the significantly enriched regions

10. For Metaplots and heatmaps, computeMatrix was used in combination with plotProfile/plotHeatmap as per instructions from deepTools web page




