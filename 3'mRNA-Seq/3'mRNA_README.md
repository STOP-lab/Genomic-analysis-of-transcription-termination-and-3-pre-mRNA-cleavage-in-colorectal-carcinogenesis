# 3'mRNA-Seq Analysis folder structure

mkdir 3Prime-Seq/FastQ 3Prime-Seq/QC 3Prime-Seq/MultiQC 3Prime-Seq/Trimmed 3Prime-Seq/Aligned 3Prime-Seq/Replicates 3Prime-Seq/Trimmed/PolyA-T/
mkdir 3Prime-Seq/QC/FastQC 3Prime-Seq/QC/TrimmedQC 3Prime-Seq/QC/Aligned 3Prime-Seq/QC/TrimmedQC/PolyA-T
mkdir 3Prime-Seq/MultiQC/FastMQC 3Prime-Seq/MultiQC/TrimmedMQC 3Prime-Seq/MultiQC/AlignedM 3Prime-Seq/MultiQC/TrimmedMQC/PolyA-T
mkdir 3Prime-Seq/Aligned/ReadCounts/
mkdir 3Prime-Seq/Aligned/Merged_ReadCounts/ 3Prime-Seq/Aligned/Stranded_BAMs/ 3Prime-Seq/Aligned/Stranded_bedGraphs/ 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs 3Prime-Seq/Aligned/Stranded_bigWigs  

1. DeMultiplexing and FastQ format conversion

# Raw .cbcl files to fastq conversion and demultiplexing  was done using custom created SampleSheet
# --no-lane-splitting make sures that the reads of one sample from different lanes don't produce multiple fastq files but would end up in one FastQ file

bcl2fastq -R 231030_A01680_0155_AHMN5MDRX3/ -o 231030_A01680_0155_AHMN5MDRX3/3Prime-Seq/FastQ/ --no-lane-splitting -p 26

2. Check the quality of reads with FastQC

fastqc 3Prime-Seq/FastQ/*.gz -o 3Prime-Seq/QC/FastQC/

multiqc  3Prime-Seq/QC/FastQC/ -n RawData_multiQC -o  3Prime-Seq/MultiQC/FastMQC/

3. Removal of adapters and PolyA/T_tails using TrimGalore => 3Prime-Seq_processing.py

multiqc 3Prime-Seq/QC/TrimmedQC/ -n TrimmedData_multiQC -o  3Prime-Seq/MultiQC/TrimmedMQC/

multiqc 3Prime-Seq/QC/TrimmedQC/PolyA-T -n PolyA-T-Trimmed_multiQC -o  3Prime-Seq/MultiQC/TrimmedMQC/PolyA-T

4. Align 3'mRNA reads to the hg38 reference genome using STAR => 3Prime-Seq_processing.py

mv 3Prime-Seq/Aligned/STARReadsPerGene.out.tab 3Prime-Seq/Aligned/ReadCounts/

mv 3Prime-Seq/Aligned/STARLog.final.out 3Prime-Seq/QC/Aligned

multiqc 3Prime-Seq/QC/Aligned/ -n 3Prime-Seq/Aligned_multiQC -o  3Prime-Seq/MultiQC/AlignedM/

5. Create strand-specfic BAMs and genome coverage files(bedGraps) => 3Prime-Seq_processing.py
## -f => includes, -F => excludes;  16 => reads mapping to reverse strand
# Forward Strand reads - 5' to 3' direction; use flag -F 16 
# Reverse Strand reads - 3' to 5' direction; use flag -f 16

## Move Normalised_bedgraphs to its respective folder
mv 3Prime-Seq/Aligned/Stranded_bedGraphs/*.bedGraph 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/

##Sort normalised bedgraphs and generate bigWigs
for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do sort -k1,1 -k2,2n "$f" -o "$f"; done

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do /home/micgdu/kentutils/bedGraphToBigWig "$f" hg38_chromsizes.genome "$f.bw"; done

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph.bw/.bw/g)"; done

mv 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw 3Prime-Seq/Aligned/Stranded_bigWigs/

6. Filtering Internal Primming Events => InternalPrimming_Mask.sh; InternalPrimming_Filtering.R

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do sort -k1,1 -k2,2n "$f" -o "$f"; done

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bedGraph; do /home/micgdu/kentutils/bedGraphToBigWig "$f" hg38_chromsizes.genome "$f.bw"; done

for f in 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph.bw/.bw/g)"; done

mv 3Prime-Seq/Aligned/Stranded_bedGraphs/Normalised_bedGraphs/*.bw 3Prime-Seq/Aligned/Stranded_bigWigs/

7. Polyadenylation site (PAS) Detection and Quantification => PAS_Calling.R; PAS_Quantification.R

######################## Renaming and relocating files after PAS Quantification ########################
for f in IPF/polyAsites/PAS/PA_counts/*.bedGraph; do mv "$f" "$(echo "$f" | sed s/_Fwd.bam_raw_PAS_Fwd.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*.bedGraph; do mv "$f" "$(echo "$f" | sed s/_Rev.bam_raw_PAS_Rev.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*FwdApa.RData; do mv "$f" "$(echo "$f" | sed s/_Fwd.bam_raw_PAS_Fwd.RData_/_raw_PAS_/)"; done
for f in IPF/polyAsites/PAS/PA_counts/*RevApa.RData; do mv "$f" "$(echo "$f" | sed s/_Rev.bam_raw_PAS_Rev.RData_/_raw_PAS_/)"; done

mv IPF/polyAsites/PAS/PA_counts/*.bedGraph IPF/polyAsites/PAS/PA_counts/bedGraphs/

######### Sorting the bedgraphs and converting them to bigWigs
for f in IPF/polyAsites/PAS/PA_counts/bedGraphs/*.bedGraph ; do sort -k1,1 -k2,2n "$f" -o "$f"; done
for f in IPF/polyAsites/PAS/PA_counts/bedGraphs/* ; do /home/micgdu/kentutils/bedGraphToBigWig "$f" /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome "$f.bw" ; done
mv IPF/polyAsites/PAS/PA_counts/bedGraphs/*.bw IPF/polyAsites/PAS/PA_counts/bigWigs/
for f in IPF/polyAsites/PAS/PA_counts/bigWigs/*.bw; do mv "$f" "$(echo "$f" | sed s/.bedGraph//)"; done

8. Alternative PolyAdenylation (APA) analysis => APA.R

9. Merge replicates of each cell line => Merge_BAMs-Replicates.sh

10. For Metaplots and heatmaps, computeMatrix was used in combination with plotProfile/plotHeatmap as per instructions from deepTools web page





