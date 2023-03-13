Projet amylose ref:

https://scholar.google.com/scholar?q=amylose+content+genetic+basis+tuber&hl=fr&as_sdt=0&as_vis=1&oi=scholart

https://www.mdpi.com/1422-0067/23/9/4640


https://www.pnas.org/doi/pdf/10.1073/pnas.2014860117


https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-020-02666-z


https://www.google.com/search?q=AGPase+gene+maize+amylose&sxsrf=AJOqlzUxP6bb4pzo9pq4dtKsHOlDaAgZmA%3A1678722469184&ei=pUUPZJPmCtmEhbIP4Ly0uAU&ved=0ahUKEwjTzM7DoNn9AhVZQkEAHWAeDVcQ4dUDCA8&uact=5&oq=AGPase+gene+maize+amylose&gs_lcp=Cgxnd3Mtd2l6LXNlcnAQAzIHCCEQoAEQCjoKCAAQRxDWBBCwA0oECEEYAFCcBFjKE2C1FmgBcAF4AIABigKIAb0OkgEDMi04mAEAoAEByAECwAEB&sclient=gws-wiz-serp



https://www.frontiersin.org/articles/10.3389/fpls.2023.1131975/full



https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-7652.2004.00073.x

















## Step 01 Grab only the genic SSR coordinates from D. alata

```bash

cat alata_genome_out_annotation.tsv | awk '/Genic/ { print }' > genic.ssr.d.alata.txt

```


## Step 02 Grab only the geneic SSR coordinates from D. roduntada

```bash


cat rotundata_genome_out_annotation.tsv | awk '/Genic/ { print }' > genic.ssr.d.rodundata.txt

```


## Step 03 find genic regions that overlap with microsatellite SSR using python code


Two inputs are required: `microsatellite_SSR_file.txt` and `genome_annotation_file.gff`

For D. rodundata I prepared the GFF file with only the 20 chromosome sets (NCBI annotation include chloroplast and mitochondrial sequences also)
For D. alata, I could not since the GFF was to big. In Nota Bene I put the sequences IDs of the 20 chromosomes of interest.


Here is the python code to grab the genic SSR with their corresponding Attributes


```python

#importing required libraries
import pandas as pd

#reading SSR file
SSR = pd.read_csv("microsatellite_SSR_file.txt", sep='\t', header=None)
SSR.columns = ['Chromosome', 'Start', 'End']

#reading GFF file
GFF = pd.read_csv("genome_annotation_file.gff", sep='\t', header=None)
GFF.columns = ['Seqid', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes']

#filtering CDS regions
CDS = GFF[GFF['Type'] == 'CDS']

#merging both files based on overlapping regions
result = pd.merge(CDS, SSR, on=['Chromosome', 'Start', 'End'], how='inner')

#extracting CDS ID
result['CDS_ID'] = result['Attributes'].str.split(";", expand = True)[0].str.split("=", expand = True)[1]

#returning final result
print("CDS Position, Start, End and Corresponding CDS ID:")
print(result[['Start', 'End', 'CDS_ID']])

```



NB: D. alata file was heavy I could not prepare the GFF file with only the 20 chromosomes sequences in excel


```

CM037011.1
CM037012.1
CM037013.1
CM037014.1
CM037015.1
CM037016.1
CM037017.1
CM037018.1
CM037019.1
CM037020.1
CM037021.1
CM037022.1
CM037023.1
CM037024.1
CM037025.1
CM037026.1
CM037027.1
CM037028.1
CM037029.1
CM037030.1

```


## Step 04 I will use vlookup to get the protein sequence ID, then grab them and perform eggNOG and or KoBas-I analyses


g


















another one






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
