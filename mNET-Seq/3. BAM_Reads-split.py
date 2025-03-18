#!usr/bin/env python
import os
for srr in open("Samples_mNET.txt"):
	sra=srr.strip("\n")
	
## Strand-specific BAMs
########### Forward Strand reads - 5' to 3' direction ###########
##Extracting the alignments of first in pair if they map to the forward strand
	os.system("samtools view -bh -@ 18 -f 64 -F 16 Replicates/Merged-BAMs/"+sra+"-Reps_merged.bam > Replicates/Stranded_BAMs/"+sra+"-Reps_merged_fwd1.bam")
##Extracting the alignments of the second in pair if they map to the reverse strand
	os.system("samtools view -bh -@ 18 -f 144 Replicates/Merged-BAMs/"+sra+"-Reps_merged.bam > Replicates/Stranded_BAMs/"+sra+"-Reps_merged_fwd2.bam")

## Merge both the fwd1 and fwd2 bams
	os.system("samtools merge -f -@ 18 Replicates/Stranded_BAMs/"+sra+"-Reps_merged_fwd.bam Replicates/Stranded_BAMs/"+sra+"-Reps_merged_fwd1.bam Replicates/Stranded_BAMs/"+sra+"-Reps_merged_fwd2.bam")
	os.system("samtools index -@ 18 -b Replicates/Stranded_BAMs/"+sra+"-Reps_merged_fwd.bam")
 
########### Reverse Strand reads - 3' to 5' direction ###########
##Extracting the alignments of first in pair if they map to the reverse strand
	os.system("samtools view -bh -@ 18 -f 80 Replicates/Merged-BAMs/"+sra+"-Reps_merged.bam > Replicates/Stranded_BAMs/"+sra+"-Reps_merged_rev1.bam")
##Extracting the alignments of the second in pair if they map to the forward strand
	os.system("samtools view -bh -@ 18 -f 128 -F 16 Replicates/Merged-BAMs/"+sra+"-Reps_merged.bam > Replicates/Stranded_BAMs/"+sra+"-Reps_merged_rev2.bam")

## Merge rev1 and rev2 bams
	os.system("samtools merge -f -@ 18 Replicates/Stranded_BAMs/"+sra+"-Reps_merged_rev.bam Replicates/Stranded_BAMs/"+sra+"-Reps_merged_rev1.bam Replicates/Stranded_BAMs/"+sra+"-Reps_merged_rev2.bam")
	os.system("samtools index -@ 18 -b Replicates/Stranded_BAMs/"+sra+"-Reps_merged_rev.bam")

## Genome Coverage files (BAM to bedGraph Conversion)
	os.system("bedtools genomecov -ibam Replicates/Stranded_BAMs/"+sra+"-Reps_merged_fwd.bam -bga -split | gzip - > Replicates/Stranded_bedGraphs/"+sra+"-Reps_merged_fwd.bedgraph.gz")
	os.system("bedtools genomecov -ibam Replicates/Stranded_BAMs/"+sra+"-Reps_merged_rev.bam -bga -split | gzip - > Replicates/Stranded_bedGraphs/"+sra+"-Reps_merged_rev.bedgraph.gz")
