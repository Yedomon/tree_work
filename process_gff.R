library(dplyr)
library(stringr)

chrom_roduntata<-as.data.frame(table(d_rodundata_gff_ready_to_use$Seqid)) 
chrom_roduntata$ChrID<-paste0("chr", rep(1:20))
colnames(chrom_roduntata)<-c("Seqid", "Freq", "ChrID")

### Merge both data 


df_chrom.merge<-left_join(d_rodundata_gff_ready_to_use, chrom_roduntata, by="Seqid") 

## Extract only protein ID from the column Attributes
df_chrom.CDS<-df_chrom.merge[df_chrom.merge$Type=="CDS",]

df_chrom.CDS.ProtID<- df_chrom.CDS %>%
                      mutate(ProtID= gsub("Name=(XP_[^;]+);$", "\\1",
                                          str_extract(Attributes, "Name=(XP_[^;]+);"))) %>%
             select(ProtID) %>% unique() %>% as.data.frame()
  

write.csv(df_chrom.CDS.ProtID, "proteinID.txt", row.names = FALSE)

## clean fasta headers
#cut -f1 -d ' '  d_rodundata_protein.faa > d_rodundata_protein.clean.faa


### extracxt the protein sequnces with ther iDs in proteinID.csv

#faSomeRecords d_rodundata_protein.clean.faa proteinID.csv out.fa



  
  
  
       #mutate(ProtID= str_extract(Attributes, "Name=(XP_[^;]+);")) %>%
       # mutate(ProtID=gsub("Name=", "", ProtID)) %>%
      #select(ProtID)









