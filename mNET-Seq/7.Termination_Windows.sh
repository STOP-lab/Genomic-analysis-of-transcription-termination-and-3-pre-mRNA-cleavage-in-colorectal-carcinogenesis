## Multiplying the reverse strand bedgraphs with -1 using macs2 bdgopt function
for f in Replicates/Stranded_bedGraphs/Normalised_bedGraps/*rev_norm.bedgraph ; do macs2 bdgopt -i "$f" -m multiply -p -1 -o "$f"-m1; done;

## Renaming reverse_norm_bedgraph files
for f in Replicates/Stranded_bedGraphs/Normalised_bedGraps/*-m1; do mv "$f" "$(echo "$f" | sed s/rev_norm.bedgraph-m1/rev_norm-m1.bedgraph/)"; done;

### Termination windows identification
## Forward strand
for f in Replicates/Stranded_bedGraphs/Normalised_bedGraps/*fwd_norm.bedgraph ; do macs2 bdgbroadcall -i "$f" -o "$f"-bdg_peaks --no-trackline ; done;
## Reverse strand
for f in Replicates/Stranded_bedGraphs/Normalised_bedGraps/*rev_norm-m1.bedgraph ; do macs2 bdgbroadcall -i "$f" -o "$f"-bdg_peaks --no-trackline ; done;

## Move files to the Termiantion_Windows folder and rename accordingly
mv Replicates/Stranded_bedGraphs/Normalised_bedGraps/*-bdg_peaks Replicates/Termination_Windows/
for f in Replicates/Termination_Windows/*_peaks.bed; do mv "$f" "$(echo "$f" | sed s/_norm.bedgraph//)"; done;
for f in Replicates/Termination_Windows/*_peaks.bed; do mv "$f" "$(echo "$f" | sed s/_norm-m1.bedgraph//)"; done;
