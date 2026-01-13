* run compleasm with dothideomycetes_odb10
  
```
#!/bin/bash

assembly=$1
lineage=$2
#library_path=/project/fdwsru_fungal/Nick/databases/lineages
threads=16

eval "$(conda shell.bash hook)" #intialize shell for conda environments
conda activate /project/fdwsru_fungal/Nick/conda/envs/compleasm

compleasm run \
        --assembly_path $assembly \
        --lineage $lineage \
        --threads $threads \
        --output_dir ${assembly%.*}_compleasm_${lineage}
```

Resutls show high completion and low duplication:

CHIP77
## lineage: dothideomycetes_odb10
S:99.71%, 3775
D:0.05%, 2
F:0.03%, 1
I:0.00%, 0
M:0.21%, 8
N:3786

EMETT52
## lineage: dothideomycetes_odb10
S:99.71%, 3775
D:0.03%, 1
F:0.08%, 3
I:0.00%, 0
M:0.18%, 7
N:3786

NGKOG
## lineage: dothideomycetes_odb10
S:99.68%, 3774
D:0.03%, 1
F:0.08%, 3
I:0.00%, 0
M:0.21%, 8
N:3786

```

Summarize results

* In folder with compleasm results

```
#!/bin/bash

echo -e "Sample\tSingle\tDuplicated\tFragmented\tIncomplete\tMissing\tTotal" > compleasm_summary.tsv


for file in */*summary.txt; do

grep -v '#' $file > tmp
sample=$(dirname $file)
S=$(awk 'NR == 1 {print}' tmp | cut -f2 -d' ')
D=$(awk 'NR == 2 {print}' tmp | cut -f2 -d' ')
F=$(awk 'NR == 3 {print}' tmp | cut -f2 -d' ')
I=$(awk 'NR == 4 {print}' tmp | cut -f2 -d' ')
M=$(awk 'NR == 5 {print}' tmp | cut -f2 -d' ')
N=$(awk 'NR == 6 {print}' tmp | cut -f2 -d':')
echo -e "$sample\t$S\t$D\t$F\t$I\t$M\t$N" >> compleasm_summary.tsv
done
rm tmp
```

Use sed to format sample names
```
Sample  Single  Duplicated      Fragmented      Incomplete      Missing Total
CHIP77  3775    2       1       0       8       3786
EMETT52 3775    1       3       0       7       3786
NGKOG   3774    1       3       0       8       3786
```


## GFAStats

Run gfastats on all fasta files, saving output as gfastats/SAMPLE.gfastats.out

summarize with script

```
#!/bin/bash
dir=$1

arr=(2 3 5 7 8 10 11 13 15 16 27)
echo "ID,scaffolds,Total_scaffold_length,Scaffold_N50,Scaffold_L50,Largest_scaffold,contigs,Total_contig_length,Contig_N50,Contig_L50,Largest_contig,GC_content_%" > gfastats_summary.csv
for file in ${dir}/*gfastats.out; do
    echo -n $(basename $file _gfastats.out),
    for i in ${arr[@]}; do
        echo -n $(sed -n "${i}p" $file | cut -f2- -d":"),
    done | sed "s/,$//1"
    echo ""
done >> gfastats_summary.csv
```

```
ID,scaffolds,Total_scaffold_length,Scaffold_N50,Scaffold_L50,Largest_scaffold,contigs,Total_contig_length,Contig_N50,Contig_L50,Largest_contig,GC_content_%
CHIP77,30,36200928,1397817,11,2169000,30,36200928,1397817,11,2169000,46.45
EMETT52,27,36467072,1534764,11,2117965,27,36467072,1534764,11,2117965,46.22
NGKOG,27,37372312,1628959,10,2519349,27,37372312,1628959,10,2519349,45.59
```


## Get telomeres present

* Use find telomeres script
e.g.

```
for file in ../*fasta; do 
echo -en "${file%.*}, "; 
python /project/fdwsru_fungal/Nick/git_repos/FindTelomeres/FindTelomeres.py $file | tail -1 | cut -f2 -d '(' | sed 's/ forward//1' | sed 's/ reverse)//1'; 
done

```



```
Sample, Forward telomeres, Reverse telomeres
CHIP77, 17, 18
EMETT52, 22, 14
NGKOG, 21, 23
```






