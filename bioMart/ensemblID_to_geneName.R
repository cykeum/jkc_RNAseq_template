#### Date: August 2023
#### Title: Conversion of ensemble IDs to common gene names using biomaRT
#### Name: Cheayeong Keum
#### Data used: RNA-seq practice round


### library loading
library(biomaRt)
library(dplyr)

### data loading
data<-read.csv("./bioMart/example_RNAseq_dataset.csv")

### Extract ensemble_id from dataset
colnames(data) # Check the column names in the data
id.data<-data$ï..Name # extract only the ids from the data

### make useMart object for getBM usage
ensembl<-useMart("ensembl", dataset = "mmusculus_gene_ensembl")

### Convert ensemble_id to gene name
annot.data<-getBM(attributes = c('ensembl_gene_id','mgi_symbol', 'external_gene_name'), filters = 'ensembl_gene_id', values = id.data, mart = ensembl)

### Left join by ensembl_gene_id
annotdf.data <- merge(
  x = as.data.frame(data),
  y =  annot.data,
  by.x = 'ï..Name',
  by.y = 'ensembl_gene_id',
  all.x = T)     ## keep all data from x (left join)

### re-order by decreasing TPM
annotdf.data<-annotdf.data[order(annotdf.data$TPM,decreasing = TRUE),]

### Change the rowname to index based on new order
rownames(annotdf.data)<-c(1:nrow(annotdf.data))

### export
write.csv(annotdf.data,'./bioMart/result biomartAnnotation.csv',row.names = FALSE)
