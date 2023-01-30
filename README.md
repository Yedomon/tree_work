## Grab only the genic SSR coordinates from D. alata

```bash

cat alata_genome_out_annotation.tsv | awk '/Genic/ { print }' > genic.ssr.d.alata.txt

```


## Grab only the geneic SSR coordinates from D. roduntada

```bash


cat rotundata_genome_out_annotation.tsv | awk '/Genic/ { print }' > genic.ssr.d.rodundata.txt

```


## find genic regions that overlap with microsatellite SSR using python code


```python



import pandas as pd

# read in microsatellite SSR data
ssr_df = pd.read_csv("microsatellite_ssr.txt", sep="\t")

# read in genome annotation file
gff_df = pd.read_csv("genome_annotation.gff", sep="\t")


# filter to only include genic regions
genic_df = gff_df[gff_df["type"] == "gene"]




# iterate through each microsatellite SSR
for index, row in ssr_df.iterrows():
  # get microsatellite SSR chromosome and positions
  chrom = row["chromosome"]
  start = row["start"]
  end = row["end"]
  
  # find genic regions that overlap with microsatellite SSR
  overlap_df = genic_df[(genic_df["chromosome"] == chrom) & (genic_df["start"] <= end) & (genic_df["end"] >= start)]
  
  # check if any overlaps were found
  if overlap_df.empty:
    # no overlaps found, microsatellite SSR is non-genic
    print("Microsatellite SSR at {}:{}-{} is non-genic".format(chrom, start, end))
  else:
    # overlaps found, microsatellite SSR is genic
    print("Microsatellite SSR at {}:{}-{} is genic".format(chrom, start, end))









```


 ## find genic regions that overlap with microsatellite SSR using R Code
 
 
 ```r
 
# Read in the microsatellite data
microsatellite_data <- read.table("microsatellite_data.txt", header = TRUE)

# Read in the genome annotation data
genome_annotation <- read.table("genome_annotation.gff", header = TRUE)

# Extract the relevant columns from the microsatellite data
microsatellite_positions <- microsatellite_data[, c("Chromosome", "Start", "End")]

# Extract the relevant columns from the genome annotation data
gene_positions <- genome_annotation[, c("Chromosome", "Start", "End")]

# Merge the two data frames by Chromosome and Start position
merged_data <- merge(microsatellite_positions, gene_positions, by = c("Chromosome", "Start"))

# Filter out the rows where the End position of the microsatellite is within the gene region
genic_microsatellites <- merged_data[merged_data$End.x <= merged_data$End.y, ]

# Filter out the rows where the End position of the microsatellite is not within the gene region
non_genic_microsatellites <- merged_data[merged_data$End.x > merged_data$End.y, ]

# View the resulting data frames
print(genic_microsatellites)
print(non_genic_microsatellites)

```








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
