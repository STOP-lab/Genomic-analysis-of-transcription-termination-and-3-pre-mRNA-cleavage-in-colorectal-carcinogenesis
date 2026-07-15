## Replace the colorectal samples with Pancreatic samples when required

# Merge the replicate bedGraph files by taking into consideration that all the genomic intervals from replicates should be covered 
bedtools unionbedg -i IPF/PAS_counts/bedGraphs/1CT_R1_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/1CT_R2_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/1CT_R3_raw_PAS_FwdApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph
bedtools unionbedg -i IPF/PAS_counts/bedGraphs/1CT_R1_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/1CT_R2_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/1CT_R3_raw_PAS_RevApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

bedtools unionbedg -i IPF/PAS_counts/bedGraphs/HCT116_R1_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/HCT116_R2_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/HCT116_R3_raw_PAS_FwdApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph
bedtools unionbedg -i IPF/PAS_counts/bedGraphs/HCT116_R1_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/HCT116_R2_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/HCT116_R3_raw_PAS_RevApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

bedtools unionbedg -i IPF/PAS_counts/bedGraphs/SW480_R1_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW480_R2_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW480_R3_raw_PAS_FwdApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph
bedtools unionbedg -i IPF/PAS_counts/bedGraphs/SW480_R1_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW480_R2_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW480_R3_raw_PAS_RevApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

bedtools unionbedg -i IPF/PAS_counts/bedGraphs/SW620_R1_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW620_R2_raw_PAS_FwdApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW620_R3_raw_PAS_FwdApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph
bedtools unionbedg -i IPF/PAS_counts/bedGraphs/SW620_R1_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW620_R2_raw_PAS_RevApaCounts_rpmC.bedGraph IPF/PAS_counts/bedGraphs/SW620_R3_raw_PAS_RevApaCounts_rpmC.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph

# Average the signal across replicates 
awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph
awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/1CT_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph
awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/HCT116_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph
awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/SW480_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_FwdApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph
awk '{total=0; for(i=4; i<=NF; i++) total+=$i; print $1, $2, $3, total/(NF-3)}' IPF/PAS_counts/bedGraphs/Replicates/SW620_unionbedg_raw_PAS_RevApaCounts_rpmc.bedGraph > IPF/PAS_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_RevApaCounts_rpmC.bedGraph

# Convert to BigWig
/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/1CT_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/1CT_Average_raw_PAS_RevApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/1CT_Average_raw_PAS_RevApaCounts_rpmC.bw

/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/HCT116_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/HCT116_Average_raw_PAS_RevApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/HCT116_Average_raw_PAS_RevApaCounts_rpmC.bw

/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/SW480_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/SW480_Average_raw_PAS_RevApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/SW480_Average_raw_PAS_RevApaCounts_rpmC.bw

/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_FwdApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/SW620_Average_raw_PAS_FwdApaCounts_rpmC.bw
/home/micgdu/kentutils/bedGraphToBigWig IPF/PAS_counts/bedGraphs/Replicates/SW620_Average_raw_PAS_RevApaCounts_rpmC.bedGraph hg38_chromsizes.genome IPF/PAS_counts/bigWigs/Replicates/SW620_Average_raw_PAS_RevApaCounts_rpmC.bw


