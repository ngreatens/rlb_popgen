filtering

* for all illumina data, run fastp, wiht adjustments for paired or single end data

```
!/bin/bash

forward_reads=$1
reverse_reads=$2
id=$3

eval "$(conda shell.bash hook)" #intialize shell for conda environments
conda activate /project/fdwsru_fungal/Nick/conda/envs/fastp
#fastp 0.23.4

fastp \
        --in1 $forward_reads \
        --in2 $reverse_reads \
        --out1 ${id}_fastp_1.fastq \
        --out2 ${id}_fastp_2.fastq \
        --html $id.html \
        --json $id.json
conda deactivate
```
or 
```
#!/bin/bash



forward_reads=$1
id=$2

eval "$(conda shell.bash hook)" #intialize shell for conda environments
conda activate /project/fdwsru_fungal/Nick/conda/envs/fastp
#fastp 0.23.4

fastp \
        --in1 $forward_reads \
        --out1 ${id}_fastp_1.fastq \
        --html $id.html \
        --json $id.json
conda deactivate
```


# Pacbio data

For pacbio datasets, run nanoplot

```
#!/bin/bash

input_fastq=$1

#intialize shell for conda environments
eval "$(conda shell.bash hook)"

conda activate /home/nicholas.greatens/.conda/envs/my_NanoPlot

NanoPlot \
        --threads 16 \
        --outdir ${1%.*}_nanoplot \
        --fastq $input_fastq

conda deactivate
```
