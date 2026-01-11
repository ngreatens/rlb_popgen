First downsampled to 30x retaining longest reads with q20 and ran with hifiasm and flye. 
Next, downsampled to 60x, retaining longest reads and assembled with  hifiasm and flye. Use with flye produced the best assemblies, with more chromosomes assembled end-to-end

Downsampling with filtlong:

```
#!/bin/bash

reads=$1

export PATH=$PATH:/project/fdwsru_fungal/Nick/git_repos/Filtlong/bin/
# Filtlong v0.2.1

filtlong \
        --min_length 10000 \
        --target_bases 2160000000 \
        $reads | gzip > ${reads%.fastq.gz}_filtered_long15_60x.fastq.gz
```

Assembly with flye

```
#!/bin/bash

reads=$1

ml flye

flye --pacbio-hifi $reads -t 32 --out-dir ${reads%.*}
```

Assemblies were assessed for telomeres with FindTelomeres.py and tidk:

```
tidk search -s TTAGGG -d out -o out assembly.fasta
tidk plot out/out_telomeric_repeat_windows.tsv
```

Tidk plots showed an assembly error in CHIP77 and EMETT52, where two contigs were joined at a telomere.
Reads were mapped to contigs and visualized in IGV for CHIP77 and EMETT52:
```
#!/bin/bash

# filter mitochrondrial reads

module load \
        minimap2 \
        samtools \
        bbtools

query=$1
target=$2

minimap2 \
        -H `#use homopolymer compressed kmer, recommended for pacbio` \
        -a `#output in sam` \
        -x map-hifi \
        -o ${query%.*}_out.sam \
        -t 72 \
        $target \
        $query

samtools view -b -F 0x100 ${query%.*}_out.sam >  ${query%.*}_primary.bam
pileup.sh in=${query%.*}_primary.bam ref=$target out=${target%.*}_pileup
```

Mapping showed assembly errors in one contig for each, due to chimeric reads. Contigs were split manually with bedtools

Next all contigs were blasted against the nt database, since some based, on pileup results has irregular coverage or GC

e.g.

1. split contigs with seqkit
```
seqkit split -i assembly.fasta
```
2. blastn all contigs against nt database
```
!/bin/bash

input_fasta=$1
db=/project/fdwsru_fungal/Nick/databases/nt/nt
threads=16


### for large assemblies, run in parallel.
## e.g. split assemblies and get list of commands "./blast.sh $assembly-part1"
# feed to parallel: "parallel < blast.cmds"
# concatenate outputs

ml blast+
blastn \
        -db $db \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -outfmt "6 qseqid staxids bitscore std sscinames stitle" \
        -num_threads $threads \
        -query $input_fasta > ${input_fasta%.*}.ncbi.blastn.out
```

All contigs were consistent with C. glycines.

Genomic contigs for each genome were subset with seqtk subseq, ordered by length, and renamed for each assembly

```
seqkit sort -l -r $assembly | awk '/^>/ {print ">contig_" ++i} !/^>/ {print}' 
```










