Spades assembly of illumina reads

```
#!/bin/bash

fwd=$1
rev=$2
outname=$3

ml spades

spades.py -t 32 -o ${outname}_spades -1 $fwd -2 $rev
```

or 

```
#!/bin/bash

reads=$1
outname=$2
ml spades


spades.py -t 32 -o ${outname}_spades -s $reads
```
