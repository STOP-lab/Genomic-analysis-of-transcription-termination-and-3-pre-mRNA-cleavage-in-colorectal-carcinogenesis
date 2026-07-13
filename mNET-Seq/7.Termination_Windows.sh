## Multiplying the reverse strand bedgraphs with -1 using the macs2 bdgopt function
for f in Replicates/Stranded_bedGraphs/Normalised_bedGraps/*rev_norm.bedgraph; do macs2 bdgopt -i "$f" -m multiply -p -1 -o "Replicates/Normalised_bedGraphs/$(basename "$f" _norm.bedgraph)_norm_m1.bedgraph"; done

### Termination windows identification
## Forward strand
for f in Replicates/Normalised_bedGraps/*fwd_norm.bedgraph ; do macs2 bdgbroadcall -i "$f" -o "Termination_Windows/$(basename "$f" _norm-m1.bedgraph)-bdg_peaks-TW" --no-trackline ; done;
## Reverse strand
for f in Replicates/Normalised_bedGraps/*rev_norm-m1.bedgraph ; do macs2 bdgbroadcall -i "$f" -o "$f"-bdg_peaks --no-trackline ; done;
