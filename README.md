# tree_work

# Tree construction with whole plastome data

We will use MAFFT (`conda install -c bioconda mafft`), TrimAl (`conda install -c bioconda trimal`) and IQ-TREE (`conda install -c bioconda iqtree`), all available on conda


# Bash script

```bash



#/bin/bash

set -e


cd /NABIC/HOME/yourworkingdirectory

### Concatenate the seven sequences

cat *.fasta > cpinput.fasta




### Multiple sequence alignment

source activate mafft_env

mafft --thread 32 --auto cpinput.fasta > cpinput.fasta.mafft

source deactivate mafft_env

### Alignment trimming

source activate trimal_env

trimal -automated1 -in cpinput.fasta.mafft -out cpinput.fasta.mafft.trimal

source deactivate trimal_env

### Tree construction

mkdir tree_construction

cd tree_construction

cp ../cpinput.fasta.mafft.trimal .

source activate iqtree_env

iqtree -s cpinput.fasta.mafft.trimal -nt AUTO -bb 1000 -alrt 1000 

source deactivate iqtree_env



## Bye!



```
