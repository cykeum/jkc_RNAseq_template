### Date: Aug 16, 2023
### Title: pipeline for RNA-seq data analysis using DESeq2 package
### Name: Cheayeong Keum

############ load libraries #############
library(DESeq2)
library(dplyr)

############ load data #############
#### 'count.raw' should be the raw (unnormalized) gene count
count.raw<-read.csv("./DEseq2/example_gene_count_matrix.csv") # plug in the path of "YOUR_FILE" here

#### 'meta.raw' contains the information about the samples
meta.raw<-read.csv("./DEseq2/example_sample_features.csv") # plug in the path of "YOUR_FILE" here

#### OPTIONAL: rename the columns to make them simple
count.raw<-count.raw %>%
  rename("GeneID"="ï..GeneID") # The former string is the name that I want. The latter is the original column name
meta.raw<-meta.raw %>%
  rename("SampleName"="ï..SampleName") # The former string is the name that I want. The latter is the original column name

############ check if the sample names in count data and metadata match ###########
all(colnames(count.raw)[-1]==meta.raw$SampleName) # Make sure this returns TRUE

#### If the test above returns FALSE, ####
#### then modify the sample names in a way that they perfectly match in the two datasets ####
#### Re-run the code above for testing again ####

#### When the test above returns TRUE, ####
#### proceed to the next step ####

############ data preprocessing ###########
#### Copy the raw data into new objects
count.modified<-count.raw
meta.modified<-meta.raw

#### Set the row name of count data to be the gene name
rownames(count.modified)<-count.raw$GeneID
count.modified<-count.modified[,-1]
head(rownames(count.modified)) # Check if the row names of count.modified are the gene names

#### Set the row name of the metadata to be the sample name
rownames(meta.modified)<-meta.raw$SampleName
meta.modified<-meta.modified[,-1]
head(rownames(meta.modified)) # Check if the row names of meta.modified are the sample names

#### OPTIONAL: Trim the unwanted / redundant columns from dataset if they exist
## In the case of example metadata I gave, ##
## the column named 'Tissue' is redundant as all samples were from placenta only ##
## So I want only the columns named 'Age' and 'Condition'
interesting<-c('Age', 'Condition')
meta.modified<-meta.modified[,interesting]

#### Convert character variables into factor
meta.modified$Age<-as.factor(meta.modified$Age)
meta.modified$Condition<-as.factor(meta.modified$Condition)

############# Prepare for DESeq input #################
#### Construct DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count.modified,
                              colData = meta.modified,
                              design = ~Age+Condition) # Variables in design must match with one of the column names in metadata


#### Pre-filtering low count genes before running DESeq2 function
smallest<-2 # Recommend to use the minimal number of samples defining discrete groups
            # Ex) In this case, we have 2 E10.5 samples. 
            #     So, I set the object 'smallest' as 2
keep<-rowSums(counts(dds) >= 10) >= smallest # Keep only rows that have a count at least 10 for minimal number of samples
dds<-dds[keep,]

#### OPTIONAL: check the unique values in factor (categorical) variables 
unique(dds$Condition) # unique values in dds$Condition are 'control' and 'treated'
unique(dds$Age) # unique values in dds$Age are 'E10.5', 'E12.5', 'E14.5', and 'E16.5'

#### re-define the levels of the factor variables
dds$Condition<-relevel(dds$Condition, ref = "control") # telling which level you want use as a reference for comparison
dds$Age<-factor(dds$Age, levels = c('E10.5', 'E12.5', 'E14.5', 'E16.5')) # explicitly setting the factor levels


################## Run DESeq #####################
dds<-DESeq(dds)

################## Result tables ####################
#### You can choose any options to see your result in this section below ####
## results() extracts a result table with log2 FC, p values and adjusted p values

#### 1. Without arguments to results() function, ####
## The comparison will be the last level of the last variable in the design formula ##
## In our case, the last variable is 'Condition' and ##
## the last level in the variable is 'treated vs control' ##
res<-results(dds)
res

#### 2. With specified coefficient or contrast we want to build a results table for ####
resultsNames(dds) # first check a list of comparisons that DESeq2 generates - it returns below
                    # [1] "Intercept"                    "Age_E12.5_vs_E10.5"          
                    # [3] "Age_E14.5_vs_E10.5"           "Age_E16.5_vs_E10.5"          
                    # [5] "Condition_treated_vs_control"

## 2.1 Build a result table for E14.5 vs E10.5 using name
res2.1<-results(dds, name = "Age_E14.5_vs_E10.5")  # name argument should be one of output from resultNames(dds) above

## 2.2 Build a result table for E14.5 vs E10.5 using contrast
res2.2<-results(dds, contrast = c("Age", "E14.5", "E10.5"))

#### 3. With specified filtering condition ####
res3_pval05 <- results(dds, alpha = 0.05) # the adjust p value cufoff is set to 0.05 now
                                          # default alpha = 0.1.

################## Summarize the result ####################
#### Here, I will proceed further with the first 'res' object

#### order our results table (for example by the smallest p value)
res_ordered<-res[order(res$pvalue),]

#### summarize some basics such as number of up and down regulated genes
summary(res)

#### filtering out only the IDs of DEGs from the result tables
deg<-rownames(res_ordered) # whole differentially expressed genes (both up and down)
up<-rownames(res_ordered[res_ordered$log2FoldChange>0,]) # a list of up-regulated genes
down<-rownames(res_ordered[res_ordered$log2FoldChange<0,]) # a list of down-regulated genes

############### Extracting deseq2 transformed values for downstream analysis ##############
vsd<-varianceStabilizingTransformation(dds, blind = FALSE)
expr_norm<-assay(vsd) # variance stablized and normalized expression profile

################# Exporting results into local files ###################
#### Export transformed values - whole gene expression profile
write.csv(as.data.frame(expr_norm),
          file = "./DESeq2/deseq normalized gene expression profile.csv")

#### Export the resulting DEG table with test statistics
write.csv(as.data.frame(res_ordered),
          file = "./DESeq2/result age E14_5 vs E10_5.csv") # Change the file name based on what you have compared with

#### Export only the list of DEG names as a text file
write.table(as.data.frame(up),
            file = "./DESeq2/upregulated genes in age E14_5 vs E10_5.txt",
            row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(down),
            file = "./DESeq2/downregulated genes in age E14_5 vs E10_5.txt",
            row.names = FALSE, col.names = FALSE)







