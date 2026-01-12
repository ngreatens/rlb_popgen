With filtered biallelic variants from 04, with low coverage and diploid/contaminated sample removed
```
VCF=SNPs_biallelic.vcf.recode.vcf
```

## first get a SNP based phylogeny

```
python /project/fdwsru_fungal/Nick/git_repos/vcf2phylip/vcf2phylip.py -i $VCF --output-folder phylogeny
cd phylogeny; ml iq_tree
iqtree2 -s SNPs_biallelic.min4.phy --seqtype DNA -m GTR
```


## prepare file for plink

Prepare rename file for contigs

```
grep '>' EMETT52.fasta | tr -d '>' | cat -n | awk -v OFS='\t' '{print $2,$1}' > rename.txt
```

Put plink2 in PATH
e.g.
```
PATH="$(echo $PATH):/90daydata/fdwsru_fungal/Nick/peanut_smut_popgen/08_popgen/"
```

```
plink2 --vcf $VCF --make-pgen --out sorted_pgen --sort-vars --rename-chrs rename.txt  --chr-set 27
plink2 --pfile sorted_pgen --make-bed --out sorted_pgen
```

* Assess population structure with admixture (submit commands to scheduler)
```
for i in {2..11}; do admixture --cv sorted_pgen.bed ${i} > log${i}.out; done
awk '/CV/ {print $3,$4}' *out 
```
Results:
```
(K=2): 0.50125
(K=3): 0.41948
(K=4): 0.37996
(K=5): 0.36845
(K=6): 0.35997
(K=7): 0.37123
(K=8): 0.38826
(K=9): 0.37552
(K=10): 0.38722
(K=11): 0.41541
```

6 has best k-value 

get names of samples from .fam file produced by plink
```
cut -f2 *fam > sample_names.txt
```

Using .Q files produced by admixture to get which pop samples belong to 

```
awk '{
    max_val = -1e9;  # Initialize max_val to a very small number
    max_col = 0;     # Initialize max_col index
    for (i = 1; i <= NF; i++) {
        if ($i > max_val) {
            max_val = $i;
            max_col = i;
        }
    }
    print max_col;
}' $1 > tmp
paste samples_names.txt tmp
```












