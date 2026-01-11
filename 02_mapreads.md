Map single and paired reads to EMETT52 reference 

```
#!/bin/bash

index=$1
fwd_reads=$2
rev_reads=$3
out_prefix=$4

ml bwa_mem2
ml samtools

bwa-mem2 mem -t 32 $index $fwd_reads $rev_reads |
samtools sort -@31 -n - |
samtools fixmate -@31  -m - ${out_prefix}_namesorted.bam
samtools sort -@31  ${out_prefix}_namesorted.bam |
samtools markdup -@31 - ${out_prefix}_markdup.bam
samtools index ${out_prefix}_markdup.bam

rm ${out_prefix}_namesorted.bam
```

or 
```
!/bin/bash

index=$1
fwd_reads=$2
out_prefix=$3

ml bwa_mem2
ml samtools

bwa-mem2 mem -t 32 $index $fwd_reads |
samtools sort -@31 -n - |
samtools fixmate -@31  -m - ${out_prefix}_namesorted.bam
samtools sort -@31  ${out_prefix}_namesorted.bam |
samtools markdup -@31 - ${out_prefix}_markdup.bam
samtools index ${out_prefix}_markdup.bam

rm ${out_prefix}_namesorted.bam
```


Following mapping, read depth was assessed with mosdepth 
```
mosdepth -n -x
```

Samples ranged in coverage from .42-61x.
