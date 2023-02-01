## Clean the fasta file header by keeping only the proteinID


```

cat your.file | awk -F' ' '{print $1}' > your_output


```


## Extract the candidate protein sequences



```

seqtk subseq your.input.fasta the_header_of_interest_IDs.list > your_output.fasta

```








# Methodology


>> The identified SSRs were classified into genic and non-genic SSRs for D. alata and D. rotuntada for which genome annotations are available. The biological attributes of genomic regions harbouring genic-SSRs were investigated with eggNOGmapper (Cantalapiedra et al., 2021) and KOBAS-i (Bu et al., 2021), respectively. 



# Code

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


## Step 04 I will use vlookup to get the protein sequence ID, 

## Step 05 Then grab them using `seqtk subseq` command 

## Step 06 Perform eggNOG and or KoBas-I analyses online
