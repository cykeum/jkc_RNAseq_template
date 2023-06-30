#### Date: Jun 28, 2023
#### Title: Conversion of ensemble IDs to common gene names using biomaRT
#### Name: Cheayeong
#### Data used: RNA-seq practice round


### library loading
library(biomaRt)
library(dplyr)

### data loading
fl <- read.csv("./rawdata/FL.csv")
jz <- read.csv("./rawdata/JZ.csv")
fp <- read.csv("./rawdata/FP.csv")
age <- read.csv("./rawdata/AGE.csv")


### filter - only top 500 highly expressed genes
fl<-fl[1:500,]
jz<-jz[1:500,]
fp<-fp[1:500,]
age<-age[1:500,]


### check duplicated gene names for manual annotations
dup.fl<-fl$Gene.name[duplicated(fl$Gene.name)]
dup.jz<-jz$Gene.name[duplicated(jz$Gene.name)]
dup.fp<-fp$Gene.name[duplicated(fp$Gene.name)]
dup.age<-age$Gene.name[duplicated(age$Gene.name)]


### Extract ensemble_id from datasets
id.fl<-fl$Name
id.jz<-jz$Name
id.fp<-fp$Name
id.age<-age$Name


### make useMart object for getBM usage
ensembl<-useMart("ensembl", dataset = "mmusculus_gene_ensembl")


### Convert ensemble_id to gene name
annot.fl<-getBM(attributes = c('ensembl_gene_id','mgi_symbol', 'external_gene_name'), filters = 'ensembl_gene_id', values = id.fl, mart = ensembl)
annot.jz<-getBM(attributes = c('ensembl_gene_id','mgi_symbol', 'external_gene_name'), filters = 'ensembl_gene_id', values = id.jz, mart = ensembl)
annot.fp<-getBM(attributes = c('ensembl_gene_id','mgi_symbol', 'external_gene_name'), filters = 'ensembl_gene_id', values = id.fp, mart = ensembl)
annot.age<-getBM(attributes = c('ensembl_gene_id','mgi_symbol', 'external_gene_name'), filters = 'ensembl_gene_id', values = id.age, mart = ensembl)


### Left join by ensembl_gene_id
annotdf.fl <- merge(
  x = as.data.frame(fl),
  y =  annot.fl,
  by.x = 'Name',
  by.y = 'ensembl_gene_id',
  all.x = T)     ## keep all data from x (left join)

annotdf.jz <- merge(
  x = as.data.frame(jz),
  y =  annot.jz,
  by.x = 'Name',
  by.y = 'ensembl_gene_id',
  all.x = T)

annotdf.fp <- merge(
  x = as.data.frame(fp),
  y =  annot.fp,
  by.x = 'Name',
  by.y = 'ensembl_gene_id',
  all.x = T)

annotdf.age <- merge(
  x = as.data.frame(age),
  y =  annot.age,
  by.x = 'Name',
  by.y = 'ensembl_gene_id',
  all.x = T)

### Change the name of the columns
resdf.fl<-annotdf.fl %>% 
  rename('MGI_symbol_manual'='MGI.symbol',
         'Gene_name_manual'='Gene.name',
         'MGI_symbol_bm'='mgi_symbol',
         'Gene_name_bm'='external_gene_name')

resdf.jz<-annotdf.jz %>% 
  rename('MGI_symbol_manual'='MGI.symbol',
         'Gene_name_manual'='Gene.name',
         'MGI_symbol_bm'='mgi_symbol',
         'Gene_name_bm'='external_gene_name')

resdf.fp<-annotdf.fp %>% 
  rename('MGI_symbol_manual'='MGI.symbol',
         'Gene_name_manual'='Gene.name',
         'MGI_symbol_bm'='mgi_symbol',
         'Gene_name_bm'='external_gene_name')

resdf.age<-annotdf.age %>% 
  rename('MGI_symbol_manual'='MGI.symbol',
         'Gene_name_manual'='Gene.name',
         'MGI_symbol_bm'='mgi_symbol',
         'Gene_name_bm'='external_gene_name')

### re-order by decreasing TPM
resdf.fl<-resdf.fl[order(resdf.fl$TPM,decreasing = TRUE),]
resdf.jz<-resdf.jz[order(resdf.jz$TPM,decreasing = TRUE),]
resdf.fp<-resdf.fp[order(resdf.fp$TPM,decreasing = TRUE),]
resdf.age<-resdf.age[order(resdf.age$TPM,decreasing = TRUE),]

### Change the rowname to index based on new order
rownames(resdf.fl)<-c(1:nrow(resdf.fl))
rownames(resdf.jz)<-c(1:nrow(resdf.jz))
rownames(resdf.fp)<-c(1:nrow(resdf.fp))
rownames(resdf.age)<-c(1:nrow(resdf.age))

### mismatch between manual annotation and biomart annotation
mismatch.fl<-resdf.fl[!(resdf.fl$MGI_symbol_manual%in%resdf.fl$MGI_symbol_bm),]
mismatch.jz<-resdf.jz[!(resdf.jz$MGI_symbol_manual%in%resdf.jz$MGI_symbol_bm),]
mismatch.fp<-resdf.fp[!(resdf.fp$MGI_symbol_manual%in%resdf.fp$MGI_symbol_bm),]
mismatch.age<-resdf.age[!(resdf.age$MGI_symbol_manual%in%resdf.age$MGI_symbol_bm),]


### mgi and gene names match for both manual annot and bm annot
# mgi_gene_mismatch.f1<-resdf.fl[!(resdf.fl$MGI_symbol_bm%in%resdf.fl$Gene_name_bm),]
# mgi_gene_mismatch.jz<-resdf.fl[!(resdf.jz$MGI_symbol_bm%in%resdf.jz$Gene_name_bm),]
# mgi_gene_mismatch.fp<-resdf.fl[!(resdf.fp$MGI_symbol_bm%in%resdf.fp$Gene_name_bm),]
# mgi_gene_mismatch.age<-resdf.fl[!(resdf.age$MGI_symbol_bm%in%resdf.age$Gene_name_bm),]


### export
write.csv(resdf.fl,'./res/biomartAnnotated_FL.csv',row.names = FALSE)
write.csv(resdf.jz,'./res/biomartAnnotated_JZ.csv',row.names = FALSE)
write.csv(resdf.fp,'./res/biomartAnnotated_FP.csv',row.names = FALSE)
write.csv(resdf.age,'./res/biomartAnnotated_AGE.csv',row.names = FALSE)


###############################################################################################
#### listing functions
# listEnsembl()
# listAttributes(ensembl)
# listFilters(ensembl)
