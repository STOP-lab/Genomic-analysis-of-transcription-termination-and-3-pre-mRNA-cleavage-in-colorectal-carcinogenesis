## Merge the replicate BAMs
samtools merge -hb Replicates/Merged-BAMs/1CT_PCF11_Reps-Merged.bam Aligned/Split_BAM/1CT_PCF11_rep1-Human.bam Aligned/Split_BAM/1CT_PCF11_rep2-Human.bam Aligned/Split_BAM/1CT_PCF11_rep3-Human.bam --threads 12
samtools merge -hb Replicates/Merged-BAMs/1CT_input_Reps-Merged.bam Aligned/Split_BAM/1CT_input_rep1-Human.bam Aligned/Split_BAM/1CT_input_rep2-Human.bam Aligned/Split_BAM/1CT_input_rep3-Human.bam --threads 12

samtools merge -hb Replicates/Merged-BAMs/HCT116_PCF11_Reps-Merged.bam Aligned/Split_BAM/HCT116_PCF11_rep1-Human.bam Aligned/Split_BAM/HCT116_PCF11_rep2-Human.bam --threads 12
samtools merge -hb Replicates/Merged-BAMs/HCT116_input_Reps-Merged.bam Aligned/Split_BAM/HCT116_input_rep1-Human.bam Aligned/Split_BAM/HCT116_input_rep2-Human.bam --threads 12

samtools merge -hb Replicates/Merged-BAMs/SW480_PCF11_Reps-Merged.bam Aligned/Split_BAM/SW480_PCF11_rep1-Human.bam Aligned/Split_BAM/SW480_PCF11_rep2-Human.bam --threads 12
samtools merge -hb Replicates/Merged-BAMs/SW480_input_Reps-Merged.bam Aligned/Split_BAM/SW480_input_rep1-Human.bam Aligned/Split_BAM/SW480_input_rep2-Human.bam --threads 12

samtools merge -hb Replicates/Merged-BAMs/SW620_PCF11_Reps-Merged.bam Aligned/Split_BAM/SW620_PCF11_rep1-Human.bam Aligned/Split_BAM/SW620_PCF11_rep2-Human.bam Aligned/Split_BAM/SW620_PCF11_rep3-Human.bam --threads 12
samtools merge -hb Replicates/Merged-BAMs/SW620_input_Reps-Merged.bam Aligned/Split_BAM/SW620_input_rep1-Human.bam Aligned/Split_BAM/SW620_input_rep2-Human.bam Aligned/Split_BAM/SW620_input_rep3-Human.bam --threads 12

## Index the merged BAMs and convert them to bigWigs
for f in Replicates/Merged-BAMs/*.bam ; do samtools index  -@ 12 -b "$f"; done;

for f in Replicates/Merged-BAMs/*.bam ; do bamCoverage -b "$f" --binSize 50 --normalizeUsing CPM -p 18 -o "$f".bw

## Peak Calling with MACS2
macs2 callpeak -t Replicates/Merged-BAMs/1CT_PCF11_Reps-Merged.bam -c Replicates/Merged-BAMs/HCEC_1CT_input_Reps-Merged.bam -f BAMPE -g hs --outdir  Replicates/Peaks/Merged-BAMs/ -n 1CT_PCF11_Reps-MergedBAM --nomodel --broad -q 0.01 --broad-cutoff 0.01 2> Replicates/Peaks/Merged-BAMs/1CT_PCF11_Reps-MergedBAM_macs2.log

macs2 callpeak -t Replicates/Merged-BAMs/HCT116_PCF11_Reps-Merged.bam -c Replicates/Merged-BAMs/HCT116_input_Reps-Merged.bam -f BAMPE -g hs --outdir  Replicates/Peaks/Merged-BAMs/ -n HCT116_PCF11_Reps-MergedBAM --nomodel --broad -q 0.01 --broad-cutoff 0.01 2> Replicates/Peaks/Merged-BAMs/HCT116_PCF11_Reps-MergedBAM_macs2.log

macs2 callpeak -t Replicates/Merged-BAMs/SW480_PCF11_Reps-Merged.bam -c Replicates/Merged-BAMs/SW480_input_Reps-Merged.bam -f BAMPE -g hs --outdir  Replicates/Peaks/Merged-BAMs/ -n SW480_PCF11_Reps-MergedBAM --nomodel --broad -q 0.01 --broad-cutoff 0.01 2> Replicates/Peaks/Merged-BAMs/SW480_PCF11_Reps-MergedBAM_macs2.log

macs2 callpeak -t Replicates/Merged-BAMs/SW620_PCF11_Reps-Merged.bam -c Replicates/Merged-BAMs/SW620_input_Reps-Merged.bam -f BAMPE -g hs --outdir  Replicates/Peaks/Merged-BAMs/ -n SW620_PCF11_Reps-MergedBAM --nomodel --broad -q 0.01 --broad-cutoff 0.01 2> Replicates/Peaks/Merged-BAMs/SW620_PCF11_Reps-MergedBAM_macs2.log
