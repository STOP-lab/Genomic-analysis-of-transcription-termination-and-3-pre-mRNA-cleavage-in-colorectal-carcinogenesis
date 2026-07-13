### Merge the BAMs of replicates
samtools merge -@ 18 -c Replicates/Merged_BAMs/SW480_DMSO.bam Aligned/SW480_DMSO_rep1_markdup.bam Aligned/SW480_DMSO_rep2_markdup.bam
samtools merge -@ 18 -c Replicates/Merged_BAMs/SW480_JTE.bam Aligned/SW480_JTE_rep1_markdup.bam Aligned/SW480_JTE_rep2_markdup.bam

samtools merge -@ 18 -c Replicates/Merged_BAMs/SW620_DMSO.bam Aligned/SW620_DMSO_rep1_markdup.bam Aligned/SW620_DMSO_rep2_markdup.bam
samtools merge -@ 18 -c Replicates/Merged_BAMs/SW620_JTE.bam Aligned/SW620_JTE_rep1_markdup.bam Aligned/SW620_JTE_rep2_markdup.bam

### Index the merged BAMs
for f in Replicates/Merged_BAMs/*.bam; do samtools index -@ 12 -b "$f"; done;
