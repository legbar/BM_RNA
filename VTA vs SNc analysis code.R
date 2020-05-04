## Now analysing VTA vs SNc neurons (hence SV)
library(DESeq2)
library(tidyverse)
library(tximeta)

setwd("Downloads/billy/")

ctsSV <- as.matrix(read.delim("Human_LCM_SNc_VTA_.txt", sep = "", row.names = "Gene_ID")) #count matrix

MetaSV <- file.path("Meta_LCM_hum_VTA_SNC.csv") # Metadata file
coldataSV <- read.csv(MetaSV, row.names = 1)
coldataSV <- coldataSV[,c("Cell_Type", "tissue")] #[] is for where the data is -- (_rows, _columns)

head(ctsSV,2)

coldataSV

rownames(coldataSV) <- sub("","", rownames(coldataSV)) 
all(rownames(coldataSV) %in% colnames(ctsSV))

all(rownames(coldataSV) == colnames(ctsSV)) # metadata names of runs == run names at header of count matrix

#ctsSV <- ctsSV[, rownames(coldataSV)] -- ignore
#all(rownames(coldataSV) == colnames(ctsSV)) -- ignore


ddsSV <- DESeqDataSetFromMatrix(countData = ctsSV, colData = coldataSV, design = ~tissue)
ddsSV

#featureData <- data.frame(tissue = rownames(ctsSV)) -- ignore 
#mcols(ddsSV) <- DataFrame(mcols(ddsSV), featureData) -- ignore
#mcols(ddsSV) -- ignore
##Addition of a feature to the metadata???? Look at the website again -- ignore


ddsSV$tissue <- relevel(ddsSV$tissue, ref = "snc" )

ddsSVDESeq<- DESeq(ddsSV)
resSV <- results(ddsSVDESeq)

resSV

resultsNames(ddsSVDESeq)
summary(resSV)

sum(resSV$padj <0.1, na.rm = TRUE)
#28 DEGs between VTA and SNc


plotCounts(ddsSV, gene = "ENSG00000112561", intgroup = "tissue") #TFEB between different tissues
plotCounts(ddsSV, gene = "ENSG00000068323", intgroup = "tissue") #TFE3 between different tissues
#No difference between SNc and VTA with TFEB because not present


#log shrinkage -- Useful for visualisation and ranking of genes. To shrink - pass dds through lfcshrink.
resultsNames(ddsSVDESeq)
resSVLFC <- lfcShrink(ddsSVDESeq, coef = "tissue_vta_vs_snc", type = "apeglm")

## Order the results table by p-value
resOrdered <- resSV[order(resSV$pvalue),]

#Here plotting gene expression changes in vta wrt snc

plotMA(resSV, ylim = c(-30,30)) 
##  plotMA(resSVLFC) -- Can use if want to visualise a little easier perhaps -- log fold


##Can interactively select the genes that you want to observe changes in using the following function then interacvtively select them
idx <- identify(resSV$baseMean, resSV$log2FoldChange)
rownames(resSV)[idx]
#Upregulated in other tissues between snc and vta/MN

#Check gene selection in graph for differential expression between samples.
plotCounts(ddsSV, gene = "ENSG00000207484", intgroup = "tissue") #Note these are DEGs between MNs and VTA - ie. what makes snc differ from both


##Export Significant (p.adj<0.1) genes into csv file
resSig <- subset(resOrdered, padj < 0.1)
resSig
write.csv(as.data.frame(resSig), file = "SNc_vs_VTA_res_sig.csv")
