# tree_work

# Tree construction with whole plastome data

We will use MAFFT, TrimAl and IQ-TREE, all available on conda


# Bash script

```bash



#/bin/bash

set -e


cd /NABIC/HOME/yedomon1/plastomics/20.phylo/set7

for i in *.faa


do


### Multiple sequence alignment

source activate mafft_env

mafft --thread 32 --auto $i > $i.mafft

source deactivate mafft_env

### Alignment trimming

source activate trimal_env

trimal -automated1 -in $i.mafft -out $i.mafft.trimal

source deactivate trimal_env


done

##---Make the matrix

cat *.mafft.trimal | awk -v RS=">" -v FS="\n" -v OFS="\n" '{for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > combined.awk.fasta

mkdir tree_construction

cd tree_construction

cp ../combined.awk.fasta .


source activate iqtree_env


iqtree -s combined.awk.fasta -nt AUTO -bb 1000 -alrt 1000 

source deactivate iqtree_env



## Bye!



bash run_phylo.sh &> log.phylo &






```
