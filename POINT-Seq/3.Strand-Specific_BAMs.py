#!usr/bin/env python
import os
for srr in open("Samples_POINT.txt"):
	sra=srr.strip("\n")

### Generating strand-specific BAMs (NEBNext Ultra II Directional RNA libraries are reverse-stranded [UTP/fr-firststrand])

### Forward transcript strand (+)
# Reads originating from transcripts on the + genomic strand are represented by: 
#- f 83: R1 mapped to the reverse strand and - f 163: R2 mapped to the forward strand

## Extracting properly paired R1 alignments mapped to the reverse strand
	os.system("samtools view -bh -@ 18 -f 83 Aligned/"+sra+"_markdup.bam > Stranded_BAMs/"+sra+"_fwd1.bam")

## Extracting properly paired R2 alignments mapped to the forward strand
	os.system("samtools view -bh -@ 18 -f 163 Aligned/"+sra+"_markdup.bam > Stranded_BAMs/"+sra+"_fwd2.bam")

## Merge both the fwd1 and fwd2 bams and index the merged BAM
	os.system("samtools merge -f -@ 18 Stranded_BAMs/"+sra+"_fwd.bam Stranded_BAMs/"+sra+"_fwd1.bam Stranded_BAMs/"+sra+"_fwd2.bam")
 	os.system("samtools index -@ 18 -b Stranded_BAMs/"+sra+"_fwd.bam")
 
### Reverse transcript strand (-)
# Reads originating from transcripts on the - genomic strand are represented by:
#- f 99: R1 mapped to the forward strand and - f 147: R2 mapped to the reverse strand (-f => includes, -F => excludes)

## Extracting properly paired R1 alignments mapped to the forward strand
	os.system("samtools view -bh -@ 18 -f 99 Aligned/"+sra+"_markdup.bam > Stranded_BAMs/"+sra+"_rev1.bam")

## Extracting properly paired R1 alignments mapped to the reverse strand
	os.system("samtools view -bh -@ 18 -f 147 Aligned/"+sra+"_markdup.bam > Stranded_BAMs/"+sra+"_rev2.bam")

## Merge rev1 and rev2 bams
	os.system("samtools merge -f -@ 18 Stranded_BAMs/"+sra+"_rev.bam Stranded_BAMs/"+sra+"_rev1.bam Stranded_BAMs/"+sra+"_rev2.bam")
	os.system("samtools index -@ 18 -b Stranded_BAMs/"+sra+"_rev.bam")

### Genome Coverage files (BAM to bedGraph Conversion)
	os.system("bedtools genomecov -ibam Stranded_BAMs/"+sra+"_fwd.bam -bga -split | gzip - > Stranded_bedGraphs/"+sra+"_fwd.bedgraph.gz")
	os.system("bedtools genomecov -ibam Stranded_BAMs/"+sra+"_rev.bam -bga -split | gzip - > Stranded_bedGraphs/"+sra+"_rev.bedgraph.gz")
