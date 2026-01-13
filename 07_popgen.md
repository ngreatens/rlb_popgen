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

Inspection of the tree shows mostly long branches, but with several very short branches typical of sexual reproduction
The following groups of samples are clones

```
KAB3_22
KAB7_22
NAM_GN_HARUN
KAB8
KAB3

KABB_2
KAB_DH_11
KAB_MR_22

NAM_5
NAM_6

KABB_1
NAL_1

NAMA_3
NAMA1

MPO_P24
MPO_P25

MPO_P12
MPO_P14

Pg_43
Pg_42

Pg_36_3
PG_36

AH_T_K1
AHTK

CHO_SS_P28
CHO_55_P28
CHO_SS_P18

Pg_34
Pg_31
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

## Assess population structure with admixture (submit commands to scheduler)
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

6 has best k-value, with 4-9 close. Four consistent with phylogeny.

get names of samples from .fam file produced by plink
```
cut -f2 *fam > sample_names.txt
```

to visualize easily with tree, make tables for each of the .Q outputs:
```
for file in *Q; do  i=$(echo $file | cut -f2 -d '.'); paste sample_names.txt $file | sed 's/ /\t/g' > admixture_pop${i}.tsv; done
```



```

## Visualize pupulation structure on tree

``` {r}
# SET WD
setwd('/90daydata/fdwsru_fungal/Nick/rlb_popgen/09_tree')


library(ggtree)


#### Format data for tree

# read in table of mat types
mat_types <- read.table('mat_types.tsv', header = TRUE)
row.names(mat_types) <- mat_types$Sample
mat_types$Sample <- NULL

# table from admixture
# prepare by pasting col 2 of .fam file from plink with admixture results and sed replacing space with tab

popdata='admixture_pop4.tsv'

pop_data <- read.table(popdata)
pop_data <- as.data.frame(pop_data, row.names = pop_data$V1)
pop_data$V1 <- NULL

#rename columns 
colnames <- c()
for (i in 1:ncol(pop_data)) {colnames <- append(colnames, i)}
colnames(pop_data) <- colnames

 

num_samples <- length(row_number(pop_data))
num_pops <- ncol(pop_data)
num_rows <- num_samples * num_pops
pop_table <- data.frame(matrix(ncol = 3, nrow = num_rows))
colnames(pop_table) <- c("Samples", "Pop", "p")
pops <- colnames

counter_1 <- 1
for (i in 1:num_samples){
  for (j in 1:num_pops){
    pop_table[counter_1, 1] <- rownames(pop_data)[i]
    pop_table[counter_1, 2] <- pops[j]
    pop_table[counter_1, 3] <- pop_data[i,j]
    counter_1 <- counter_1 + 1
  }
}

pop_table$Pop <- as.factor(pop_table$Pop)

# merge tables 
df <- merge(mat_types, pop_data, by = "row.names")
df$MAT <- as.factor(df$MAT)


#read in tree
tree <- read.tree("tree.nwk")

#midpoint root tree 
library(phytools)
library(ggplot2)
library(ggstance)
tree <- midpoint.root(tree)

textsize=1.5

p <- ggtree(tree)  %<+% df +
     geom_tiplab(
      aes(color = MAT), #format excel with a Color column
      parse = TRUE,
      align = TRUE,
      size = textsize, key_glyph = "rect"
    ) +  xlim_expand(.3, 'Tree') +
  geom_facet(
    panel = 'pops', 
    data = pop_table, 
	  geom = geom_barh, 
  	mapping = aes(x = p, fill = Pop), 
    stat='identity', width = .8) +
	  theme_tree2(
	    strip.background = element_blank(), 
	    strip.text.x = element_blank()
	    ) + xlim_expand(4, 'pops')
```



```
mkdir 6_pop

for i in {1..6}; do 
	mkdir 6_pop/list_${i}
	awk -v var="$i" '($2 == var) {print $1}' tmp2 > 6_pop/list_${i}/list
done

for i in {1..6}; do 
	bcftools view -S 6_pop/list_${i}/list SNPs_biallelic.vcf.recode.vcf > 6_pop/list_${i}/pop_${i}.vcf
done

for i in {1..6}; do 
	cd 6_pop/list_${i}
	plink2 --vcf pop_${i}.vcf --make-pgen --out sorted_pgen --sort-vars --rename-chrs /90daydata/fdwsru_fungal/Nick/rlb_popgen/07_popgen/rename.txt  --chr-set 27 --mind .5 --geno 0.1 --vcf-half-call missing --maf .2
	plink2         -pfile sorted_pgen         --ld-window 100         --ld-window-kb 1000         --r2-phased cols=+dprimeabs         --ld-window-r2 0
	cd ../../
done
```

* Run R script for each plink2.vcor

```
df <- read.table('plink2.vcor', header = FALSE)
BINSIZE=1000
df$dist <- df$V5 - df$V2
df$bin <- round(df$dist/BINSIZE, 0)

library(plyr)
df2 <- ddply(df, .(bin), summarise,
      meanr2 = mean(V7))
df3 <- ddply(df, .(bin), summarise,
      meandprime = mean(V8))

plot(df2$bin*BINSIZE/1000, df2$meanr2, xlab="Physical distance (kbp)", ylab="r-squared", main="r-squared decay rate", ylim = c(0, 1))
plot(df2$bin*BINSIZE/1000, df3$meandprime, xlab="Physical distance (kbp)", ylab="D'", main="D' decay rate", ylim = c(0, 1))


png(filename = "r2_vs_dist.png", width = 400, height = 300, units = "px")
plot(df2$bin*BINSIZE/1000, df2$meanr2, xlab="Physical distance (kbp)", ylab="r-squared", main="r-squared decay rate", ylim = c(0, 1))
dev.off()

png(filename = "Dprime_vs_dist.png", width = 400, height = 300, units = "px")
plot(df2$bin*BINSIZE/1000, df3$meandprime, xlab="Physical distance (kbp)", ylab="D'", main="D' decay rate", ylim = c(0, 1))
dev.off()
```

Everything summarized with images in a powerpoint file in rlb_phylogenomic 01_13_26...










