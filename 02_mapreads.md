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

Add RGs
```
#!/bin/bash

module load \
        java/17 \
        picard/3.0.0


input=$1
id=$2
output=${1%.*}_new.bam

#library
RGLB=${id}_1

#platform
RGPL=ILLUMINAHISEQ2

#sample_name
RGSM=${id}

#ID
RGID=${id}

#platform unit e.g. run barcode
RGPU=${RGID}_001

java -Xmx16G -jar /software/el9/apps/picard/3.0.0/picard.jar \
\
AddOrReplaceReadGroups \
        --INPUT $input\
        --OUTPUT $output \
        --RGLB $RGLB \
        --RGPL $RGPL \
        --RGSM $RGSM \
        --RGID $RGID \
        --RGPU $RGPU


ml samtools

samtools index ${1%.*}_new.bam
```
