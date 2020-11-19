# Title     : homework 1
# Objective : What genes and pathways are upregulated in untreated samples?
# Created by: ALEX
# Created on: 17.11.2020

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("sva", quietly = TRUE)) BiocManager::install("sva")
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")

library(doParallel)
#Find out how many cores are available (if you don't already know)
cores <- detectCores()
#Create cluster with desired number of cores, leave one open for the machine
#core processes
cl <- makeCluster(cores[1]-1)
#Register cluster
registerDoParallel(cl)


library(GEOquery)
library(Biobase)
library(ggplot2)
library(reshape2)
library(limma)
library(MASS)
library(fgsea)
library(pheatmap)
library(sva)
library(ggrepel)
library(dplyr)

GSE53986 <- getGEO("GSE53986", AnnotGPL = TRUE)[[1]]


pData(GSE53986)$rep <- gsub(".*, (\\d)$", "rep \\1", pData(GSE53986)$title)
pData(GSE53986) <- pData(GSE53986)[, c("cell type:ch1", "treatment:ch1", "rep")]
colnames(pData(GSE53986)) <- c("Cell", "Treatment", "Replicates")


#
pData(GSE53986)$Untreated <- as.factor(c("no", "yes")[grepl("Untreated", pData(GSE53986)$Treatment) + 1])
fData(GSE53986) <- fData(GSE53986)[, c("ID", "Gene symbol", "Gene ID")]

exprs(GSE53986) <- normalizeBetweenArrays(log2(exprs(GSE53986) + 1), method="quantile")

GSE53986 <- GSE53986[!grepl("///", fData(GSE53986)$'Gene symbol'), ]
GSE53986 <- GSE53986[fData(GSE53986)$'Gene symbol' != "", ]

fData(GSE53986)$mean_expression <- apply(exprs(GSE53986), 1, mean)
GSE53986 <- GSE53986[order(fData(GSE53986)$mean_expression, decreasing = TRUE), ]
GSE53986 <- GSE53986[!duplicated(fData(GSE53986)$'Gene symbol'), ]
GSE53986 <- GSE53986[seq_len(12000), ]

dim(GSE53986)


# PCA
pcas <- prcomp(t(exprs(GSE53986)), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(GSE53986))

ggplot(plotData, aes(x=PC1, y=PC2, color=Untreated, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)

GSE53986.design <- model.matrix(~1 + Untreated, data=pData(GSE53986))

fit <- lmFit(GSE53986, GSE53986.design)

fit2 <- contrasts.fit(fit, makeContrasts(Untreatedyes, levels=GSE53986.design))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method = "BH", number=Inf, sort.by = "P")

ggplot(de, aes(x=logFC, y=-log10(adj.P.Val), color=adj.P.Val < 0.05)) +
  geom_point() + theme_bw() +
  geom_text(data = de[de$adj.P.Val < 0.01,], color = 'black', size = 2,   # and here
            aes(label=Gene.symbol))

head(de)

upRegulatedGenes <- de %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull('Gene.symbol')
length(upRegulatedGenes)

# Pathway enrichment
load("D:/R_files/keggSymbolMouse.rdata")

stats <- de$t
names(stats) <- de$Gene.symbol

fgseaResults <- fgseaMultilevel(keggSymbolMouse, stats, minSize = 15, maxSize = 500)

topPathwaysUp <- fgseaResults[ES > 0, ][head(order(pval), n=5), pathway]
topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(keggSymbolMouse[topPathways], stats, fgseaResults, gseaParam = 0.5)

names(stats)
