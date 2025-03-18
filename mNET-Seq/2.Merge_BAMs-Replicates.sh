### Merge the BAM files of replicates

## T4ph mNET-Seq samples
samtools merge -@ 14 -bh Replicates/Merged-BAMs/1CT-T4ph-Reps_merged.bam Aligned/1CT-T4ph_rep1-STARAligned.sortedByCoord.out.bam Aligned/1CT-T4ph_rep2-STARAligned.sortedByCoord.out.bam
samtools merge -@ 14 -bh Replicates/Merged-BAMs/HCT116-T4ph-Reps_merged.bam Aligned/HCT116-T4ph_rep1-STARAligned.sortedByCoord.out.bam Aligned/HCT116-T4ph_rep2-STARAligned.sortedByCoord.out.bam
samtools merge -@ 14 -bh Replicates/Merged-BAMs/SW480-T4ph-Reps_merged.bam Aligned/SW480-T4ph_rep1-STARAligned.sortedByCoord.out.bam Aligned/SW480-T4ph_rep2-STARAligned.sortedByCoord.out.bam
samtools merge -@ 14 -bh Replicates/Merged-BAMs/SW620-T4ph-Reps_merged.bam Aligned/SW620-T4ph_rep1-STARAligned.sortedByCoord.out.bam Aligned/SW620-T4ph_rep2-STARAligned.sortedByCoord.out.bam Aligned/SW620-T4ph_rep3-STARAligned.sortedByCoord.out.bam

## Total mNET-Seq samples
samtools merge -@ 14 -bh Replicates/Merged-BAMs/1CT-Total-Reps_merged.bam Aligned/1CT-Total_rep1-STARAligned.sortedByCoord.out.bam.bam Aligned/1CT-Total_rep2-STARAligned.sortedByCoord.out.bam Aligned/1CT-Total_rep3-STARAligned.sortedByCoord.out.bam
samtools merge -@ 14 -bh Replicates/Merged-BAMs/HCT116-Total-Reps_merged.bam Aligned/HCT116-Total_rep1-STARAligned.sortedByCoord.out.bam Aligned/HCT116-Total_rep2-STARAligned.sortedByCoord.out.bam
samtools merge -@ 14 -bh Replicates/Merged-BAMs/SW480-Total-Reps_merged.bam Aligned/SW480-Total_rep1-STARAligned.sortedByCoord.out.bam Aligned/SW480-Total_rep2-STARAligned.sortedByCoord.out.bam 
samtools merge -@ 14 -bh Replicates/Merged-BAMs/SW620-Total-Reps_merged.bam Aligned/SW620-Total_rep1-STARAligned.sortedByCoord.out.bam Aligned/SW620-Total_rep2-STARAligned.sortedByCoord.out.bam Aligned/SW620-Total_rep3-STARAligned.sortedByCoord.out.bam

### Index the merged BAMs
for f in Replicates/Merged-BAMs/*.bam; do samtools index -@ 12 -b "$f"; done;
