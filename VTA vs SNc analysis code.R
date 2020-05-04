## Now analysing VTA vs SNc neurons (hence SV)
library(DESeq2)
library(tximeta)

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)

getwd()

ctsSV <- as.matrix(read.delim("Human_LCM_SNc_VTA_.txt", sep = "", row.names = "Gene_ID")) #count matrix

#How many genes have >= 1 RPKM mean across samples?
sum(base::rowMeans(ctsSV) >= 1)
# >= 0.1
sum(base::rowMeans(ctsSV) >= 0.1)
# == 0 ?!?
sum(base::rowMeans(ctsSV) == 0)

#When looking at the data using RPKM, there is some serious skew.
base::summary(rowMeans(ctsSV))

#Let's plot this: Not a very clear plot initially, 
#because there are so many zero values in the data, 
#and it has not been log2 transformed
ctsSV %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample, value = count, -gene) %>%
  inner_join(coldataSV, by = c("sample" = "sample")) %>%
  ggplot(aes(count, fill = sample, colour = tissue)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw()

#remove genes <1 RPKM mean
ctsSV <- ctsSV[base::rowMeans(ctsSV) >= 1,]
base::summary(rowMeans(ctsSV))

#Still hard to visualise
ctsSV %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample, value = count, -gene) %>%
  inner_join(coldataSV, by = c("sample" = "sample")) %>%
  ggplot(aes(count, fill = sample, colour = tissue)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw()

#Let's see what log2-transformed data looks like
ctsSV_log <- log2(ctsSV) #log2(ctsSV + 1) if zero values are present in the matrix

log2(ctsSV) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample, value = count, -gene) %>%
  inner_join(coldataSV, by = c("sample" = "sample")) %>%
  ggplot(aes(count, fill = sample, colour = tissue)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw()

#Why does this matter? Mean vs variance relationship - This is where I am up to
ctsSV %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample, value = count, -gene) %>%
  group_by(gene) %>%
  mutate(mean_count = mean(count))
  

ctsSV %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  gather(key = sample, value = count, -gene) %>%
  inner_join(coldataSV, by = c("sample" = "sample")) %>%
  ggplot(aes(count, fill = sample, colour = tissue)) +
  geom_line(stat = "density", adjust = 0.75) +
  theme_bw() +
  # scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Density of Read Counts") +
  xlab("RPKM") +
  facet_grid(rows = vars(tissue)) +
  theme(legend.position = "none") +
  theme(strip.text.y = element_text(angle = 360)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#We have 65000 'genes' and 12 samples before filtering, 
#12876 after filtering out genes with less than 1 RPKM mean across samples
dim(as.matrix(ctsSV)) 

# MetaSV <- file.path("Meta_LCM_hum_VTA_SNC.csv") # Metadata file
# coldataSV <- read.csv(MetaSV, row.names = 1)

# read_delim to make a 'tidy' tibble
# Pipe the full coldataSV tibble into a 'select' function, 
# to just select the columns we want, while also renaming the capitalised names.
# Arranging the tibble by tissue
coldataSV <- read_delim("Meta_LCM_hum_VTA_SNC.csv", delim = ",") %>%
  select("sample" = Run, "cell_type" = Cell_Type, tissue) %>%
  arrange(tissue)

# coldataSV <- coldataSV[,c("Cell_Type", "tissue")] #[] is for where the data is -- (_rows, _columns)

coldataSV

head(ctsSV,2)

# rownames(coldataSV) <- sub("","", rownames(coldataSV)) 
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

#Find out what this dispersion message is about
plotDispEsts(ddsSVDESeq)


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
