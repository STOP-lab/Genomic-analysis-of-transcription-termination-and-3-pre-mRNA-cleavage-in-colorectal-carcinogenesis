# Merge the replicate bedGraph files by taking into consideration that all the genomic intervals from replicates should be covered 

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/1CT_R1_S28_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/1CT_R2_S29_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/1CT_R3_S30_raw_PAS_FwdApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/1CT_R1_S28_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/1CT_R2_S29_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/1CT_R3_S30_raw_PAS_RevApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/HCT116_R1_S31_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/HCT116_R2_S32_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/HCT116_R3_S33_raw_PAS_FwdApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/HCT116_R1_S31_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/HCT116_R2_S32_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/HCT116_R3_S33_raw_PAS_RevApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW620_R2_S25_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW620_R3_S26_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW620_R4_S27_raw_PAS_FwdApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW620_R2_S25_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW620_R3_S26_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW620_R4_S27_raw_PAS_RevApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW480_R1_S34_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW480_R2_S35_raw_PAS_FwdApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW480_R3_S36_raw_PAS_FwdApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph

bedtools unionbedg -i 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW480_R1_S34_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW480_R2_S35_raw_PAS_RevApaCounts_rpmC.bedGraph 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/SW480_R3_S36_raw_PAS_RevApaCounts_rpmC.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

# Average the signal across replicates 

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

# Convert to BigWig

/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_RevApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_RevApaCounts_rpmC.bw

/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_RevApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_RevApaCounts_rpmC.bw

/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_RevApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_RevApaCounts_rpmC.bw

/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_RevApaCounts_rpmC.bedGraph /dysk2/groupFolders/deepshika/GenomicData/hg38_chromsizes.genome 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_RevApaCounts_rpmC.bw

mv 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bedGraphs/Replicates/*.bw 3Prime-Seq/IPF/polyAsites/PAS/PA_counts/bigWigs/Replicates/

