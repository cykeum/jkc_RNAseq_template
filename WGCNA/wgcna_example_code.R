### Date: Aug 19, 2023
### Title: pipeline for weighted gene coexpression network analysis (WGCNA)
### Name: Cheayeong Keum


############ load libraries #############
library(WGCNA)
library(dplyr)
library(CorLevelPlot)

############ load data #############
#### Expression profile
expr <- read.csv("./WGCNA/example_deseq-normalized_countdata.csv")

#### Sample trait
pheno <- read.csv("./WGCNA/example_phenodata.csv")

############ my function #############
#### This is the function I made to simplify the process of network construction
network_construction<-function(expressionSet, sft_power, mergeThreshold) {
  net<-blockwiseModules(expressionSet, power = sft_power,
                        TOMType = "signed",
                        minModuleSize = 30, # relatively large min module size
                        reassignThreshold = 0,
                        mergeCutHeight = mergeThreshold, # threshold for merging of modules
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,
                        # saveTOMFileBase = "TOM",
                        verbose = 3,
                        randomSeed = 10)
  return(net)
}

############ input preparation ##############
#### Input needs to be a transposed count matrix that has been normalized by DESeq
rownames(expr)<-expr$X # change the row name as ensemble ids
expr<-expr[,-1]
input_mat <- as.data.frame(t(expr)) # transpose the count matrix

#### preprocess the phenodata into a right format and trim unnecessary ones
pheno<-pheno %>%
  select(-c('geo_id'))
rownames(pheno)<-pheno$X
pheno<-pheno[,-1]

############ Data quality check / quality control ##############
#### filtering genes with too many missing values
gsg <- goodSamplesGenes(input_mat)
summary(gsg)
gsg$allOK 

#### if gsg$allOK is true, then you can skip this step and proceed to the next
#### Otherwise, please comment out the following and run the code below
# input_mat<-input_mat[,gsg$goodGenes = TRUE]

#### Outlier detection using hierarchical clustering
htree <- hclust(dist(input_mat), method = "average") # clustering to see if there is obvious outlier

## draw a tree
par(cex = 0.6)
plot(htree) # for this example dataset, you should not see any outliers

#### OPTIONAL: If there is some obvious outliers, you can comment out the following and run the codes below ####
#### But please note that this will compensate the sample size.####
#### Ex) lets say that you want to remove the branches that are higher than height 40
# abline(h=40, col='red') # now you will see a line at the height of 2.9 on the dendrogram
# clust <- cutreeStatic(htree, cutHeight = 40, minSize = 10) # select the clusters under the line
# table(clust) # this will show outliers and samples that we want to keep
# keepSam<-(clust==1) # samples that I want to keep
# input_mat<-input_mat[keepSam,]

############# Soft power detection for network construction ################
#### soft power series generation
powers <- c(c(1:10), seq(from=12, to = 20, by=2)) # a set of thresholding powers
sft <- pickSoftThreshold(input_mat, powerVector = powers, verbose = 5) # network topology analysis

#### power visualization
par(mfrow=c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,
     cex = 0.8,
     col = "red")
abline(h=0.8, col="blue") # line corresponds to using selected an R^2 cut-off of h
plot(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity graph")
text(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     labels = powers,
     cex = 0.8,
     col = "darkgreen")

#### Based on the visualization you did above, please select the smallest value where the graphs become almost flat
#### i.e., pick a soft threshold power near the curve of the plot
#### This value will be a soft-power that we are going to use for network construction
#### Here, I set my soft-power to be 12
my_power = 12 # use the power determined from previous steps
eset <- input_mat # expression set
threshold <- 0.25 # The higher the threshold, the larger module will be generated

############# Network construction ############
#### please look at the "my function" section above if you want to adjust detailed parameters for network construction
net <- network_construction(eset, my_power, threshold)

#### show the number of genes in each module number
table(net$colors) # In this case there are 63 modules. These are sorted in decreasing order of module size

#### convert module number to a colour for each gene
mergedColors <- labels2colors(net$colors)

#### save the network
moduleLabels <- net$colors # module labels by numbers
moduleColors <- labels2colors(net$colors) # module labels by colours
MEs_label <- net$MEs # module eigengenes
geneTree = net$dendrograms[[1]] # network dendrogram
save(MEs_label, moduleLabels, moduleColors, geneTree,
     file = "./WGCNA/networkConstruction-auto.RData") # Change the file name when you save your own file

############# Identifying important modules ############
#### Define numbers of genes and samples
num_genes<- ncol(input_mat)
num_samples<-nrow(input_mat)

#### Recalculate module eigengenes with color labels
MEs_colour <- moduleEigengenes(input_mat, moduleColors)$eigengenes
MEs_colour <- orderMEs(MEs_colour) # re-order the module eigengenes

#### Create a model matrix for each trait ####
#### You need to manually do this step for all (non-numeric) traits that you are interested in ####
#### In the example, there are 2 traits: dex (design experiment) and celltype ####
dex.mod <- model.matrix(~0+pheno[['dex']]) # model matrix for trait 'dex'
colnames(dex.mod)<-levels(as.factor(pheno[['dex']]))

celltype.mod <- model.matrix(~0+pheno[['celltype']]) # model matrix for trait 'celltype'
colnames(celltype.mod)<-levels(as.factor(pheno[['celltype']]))

trait_mat <- cbind(dex.mod, celltype.mod) # create the combinatory matrix for all traits you are interested in
rownames(trait_mat)<-rownames(pheno)

#### Module to trait association
module_trait_corr<-cor(MEs_colour, trait_mat, use = "p") # p stands for pearson correlation
                                                    # this generates a table for module to trait correlation
module_trait_pval <- corPvalueStudent(module_trait_corr, num_samples) # module to trait pvalue

#### Visualize the module-trait association with a heatmap ####
#### Each cell in the heatmap contains the corresponding correlation and p-values ####
#### If a module is with 3 asterisks (***) for a certain trait, ####
####  then this module is highly likely to be a significant module for the trait ####
myheatmap <- merge(MEs_colour, trait_mat, by='row.names')
myheatmap<-myheatmap %>%
  column_to_rownames(var = "Row.names")
CorLevelPlot(myheatmap,
             x=names(myheatmap)[64:69], # select columns corresponding to traits in myheatmap
             y=names(myheatmap)[1:63], # select columns corresponding to modules in myheatmap
             main = 'Module-trait association',
             cexMain = 1.0,
             cexLabX = 0.7,
             cexLabY = 0.5,
             cexCorval = 0.5,
             rotLabX = 15,
             col = c("blue1","skyblue","white","pink","red"))

## ^ EXAMPLE: MEblue has a significant association with trait "control" (and "treated") according to the heatmap
## I will use this example (module: MEblue, and trait: control) for following steps of analysis

############# Gene relationship to trait and important modules ############
###### You can repeat this section for each of the traits that you want to look for #######
#### Pick a trait of interest - here, I chose "control"
interesting_trait <- as.data.frame(trait_mat)$control # change the dataframe column to be the one you are interested in
interesting_trait <- as.data.frame(interesting_trait)
names(interesting_trait) <- "control" # change the names of "interesting_trait" dataframe to be the trait name that you are interested in

#### Extract the names (colors) of the modules
moduleNames <- substring(names(MEs_colour),3)

#### Define a module membership as the correlation of the module eigengene and the gene expression profile ####
#### This allows the quanification of the similarity of all genes on the array to every module
geneMM <- as.data.frame(cor(input_mat, MEs_colour, use = "p")) # module membership as correlation
MMPval <- as.data.frame(corPvalueStudent(as.matrix(geneMM), num_samples)) # module membership p-values
names(geneMM)<-paste("MM", moduleNames, sep = "")
names(MMPval)<-paste("p.MM", moduleNames, sep = "")

#### Quantify associations of individual genes with a trait of interest (control) ####
####  by defining gene significance as the correlation between the gene and the trait ####
gene_trait_sig<-as.data.frame(cor(input_mat, interesting_trait, use = "p")) # Gene significance as correlation
GSPval <- as.data.frame(corPvalueStudent(as.matrix(gene_trait_sig), num_samples))
names(gene_trait_sig)<-paste("GS", names(interesting_trait), sep = "")
names(MMPval)<-paste("p.GS", names(interesting_trait), sep = "")


############# Intramodular analysis - identify genes with high GS and MM ############
###### You can repeat this section for each of the modules that you want to look for #######
#### pick a module of interest - here, I chose "blue" module
interesting_module <- "blue"
module_column <- match(interesting_module, moduleNames)
module_genes<-(moduleColors == interesting_module)

#### Visualizing module membership vs gene significance
par(mfrow=c(1,1))
verboseScatterplot(x = abs(geneMM[module_genes, module_column]),
                   y = abs(gene_trait_sig[module_genes, 1]),
                   xlab = paste("Module membership in", interesting_module, "module"),
                   ylab = paste("Gene significance for",names(interesting_trait)),
                   main = "MM vs. GS\n ",
                   cex.main = 1.0,
                   cex.lab = 1.0,
                   cex.axis = 0.9,
                   col = interesting_module
) # This graph must show higher correlation between MM vs GS for further analysis

#### If you want to look at which genes (IDs) are in the module of interest, please use the following codes
geneID_in_interestingModule <- names(input_mat)[moduleColors==interesting_module]
geneID_in_interestingModule # this returns IDs belonging to the module of interest (blue)

############# Save the summary output of network analysis results #################
###### You can repeat this section for each of the traits that you want to look for
#### Create the starting data frame - for a trait of interest
geneInfo = data.frame(geneID = names(input_mat),
                      moduleColor = moduleColors,
                      geneTraitSignificance = gene_trait_sig,
                      GSPvalue = GSPval)

#### Order modules by their significance for the trait
moduleOrder = order(-abs(cor(MEs_colour, interesting_trait, use = "p")))

#### Add module membership information in the chosen order
for (mod in 1:ncol(geneMM)) {
  oldNames = names(geneInfo)
  geneInfo = data.frame(geneInfo, geneMM[, moduleOrder[mod]],
                         MMPval[, moduleOrder[mod]]);
  names(geneInfo) = c(oldNames, paste("MM.", moduleNames[moduleOrder[mod]], sep=""),
                       paste("p.MM.", moduleNames[moduleOrder[mod]], sep=""))
}

#### Order the genes in teh geneInfo variable first by module color, then by gene-trait significance
geneOrder = order(geneInfo$moduleColor, -abs(geneInfo$GScontrol))
geneInfo <- geneInfo[geneOrder,]

#### Extract the dataframe into a spreadsheet
write.csv(geneInfo, file = "./WGCNA/geneInfo_for_Control.csv") # change the file name if you need

############# Export gene network for external visualization software called VisANT #################
#### calculate topological overlap
tom <- TOMsimilarityFromExpr(input_mat, power=my_power)

#### select a module that you want to visualize its network
mod <- "blue" # change this to be a specific module. ex) "brown"

#### Select the module genes
geneID.whole<-names(input_mat)
in_module <- (moduleColors == mod)
geneID.mod <- geneID.whole[in_module]

#### Select the corresponding topological overlap
modTom <- tom[in_module, in_module]
dimnames(modTom) <- list(geneID.mod, geneID.mod)

#### Export the network into an edge list file VisANT can read
## OPTIONAL: restrict the genes in the output to top hub genes (EX. 30 top hub genes)
n_top = 30
IMConn <- softConnectivity(input_mat[, geneID.mod])
top <- (rank(-IMConn) <= n_top)

## Export vis file
vis <- exportNetworkToVisANT(modTom,
                             file = paste("./WGCNA/VisANTInput",mod,".txt", sep=""),
                             weighted = TRUE,
                             threshold = 0)


######################################################################################
# For visualization, install visANT software (http://www.visantnet.org/visantnet.html)
# Load the file produced by the above code in VisANT.
# Then, play with various threshold and graphic layouts for better picture






