Samples with greater than 5x coverage were included in variant calling

Eight samples were excluded based on coverage. Additionally, EMETT_Q was excluded based on the presence of both MAT genes, indicating either diploidy or contamination

excluded samples
```
Pg_35
PLDT7_HARUN
NAMA_1
LW_P13
LUW_MAB_GN
KAB8_22
KAB_6
ETEP1_132
EMETT_Q
```

The remaining 125 samples were passed to freebayes for variant calling

```
#!/bin/bash

ref=$1
bamlist=$2
out=$3

module load freebayes

freebayes-parallel \
        <(fasta_generate_regions.py ${ref}.fai 100000) 144 \
        --fasta-reference ${ref} \
        --ploidy 1 \
        --bam-list $bamlist > $out
```

And filter

```
#!/bin/bash

vcf=$1

eval "$(conda shell.bash hook)" #intialize shell for conda environments
conda activate vcflib

vcffilter -f "QUAL > 20 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & AC > 0" $1 > out_filtered.vcf



ml vcftools

# separate indels
vcftools --vcf out_filtered.vcf --keep-only-indels --recode --recode-INFO-all --out output_indels-only.vcf
# separate SNPs
vcftools --vcf out_filtered.vcf --remove-indels --recode --recode-INFO-all --out output_snps-only.vcf


ml vcftools

vcftools \
        --vcf output_snps-only.vcf.recode.vcf \
        --out SNPs_biallelic.vcf \
        --min-alleles 2 \
        --max-alleles 2 \
        --max-missing .9 \
        --minDP 5 \
        --maf .05 \
        --recode
```


