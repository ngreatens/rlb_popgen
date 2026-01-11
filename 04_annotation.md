The three pacbio genomes were annotated with EDTA and Braker3

EDTA repeat annotation and masking
```
#!/bin/bash

genome=$1

cd $(dirname $genome)

eval "$(conda shell.bash hook)" #intialize shell for conda environments
conda activate /project/fdwsru_fungal/Nick/conda/envs/EDTA2

/project/fdwsru_fungal/Nick/git_repos/EDTA/EDTA.pl \
        --genome $genome \
        --species others \
        --step all \
        --threads 16 \
        --sensitive 1 \
        --overwrite 0 \
        --anno 1 \
        --force 1

#soft mask genome
perl /project/fdwsru_fungal/Nick/git_repos/EDTA/util/make_masked.pl \
        -genome $genome \
        -hardmask 0 \
        -minlen 80 \
        -rmout ${genome}.mod.EDTA.anno/${genome}.mod.EDTA.RM.out

mv ${genome}.new.masked ${genome}.masked
```
Braker 3 gene prediction

```
#!/bin/bash

masked_genome=$1
species=$2
prot_db=/project/fdwsru_fungal/Nick/databases/ortho_db/orthodb/dothideomycetes_odb11.fa

mkdir $species; cd $species

module load apptainer

apptainer exec /project/fdwsru_fungal/Nick/sifs/braker3_latest.sif braker.pl \
        --genome $masked_genome \
        --fungus \
        --species $species \
        --prot_seq $prot_db \
        --threads 32 \
        --gff3 \
        --AUGUSTUS_CONFIG_PATH /project/fdwsru_fungal/Nick/software/augustus_config \
        --GENEMARK_PATH=/project/software/el9/apps/genemark/4.71/gmes_linux_64_4
```
