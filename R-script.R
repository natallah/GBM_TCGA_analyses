library(gsean)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(WGCNA)
library(RANKS)
library(stringr)
library(data.table)

setwd("/depot/pccr/data/Matosevic/tcga_analysis/Analysis/1_gsean/GBM")

CancerProject <- "TCGA-GBM"

query <- GDCquery(project = CancerProject,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = T)

GDCdownload(query, method = "api")

invisible(capture.output(data <- GDCprepare(query)))

exprs.GBM <- assay(data)
#saveRDS(exprs.GBM, file = "exprs.GBM.rds")

exprs.GBM = readRDS("exprs.GBM.rds")


# remove duplicated gene names
exprs.GBM <- exprs.GBM[-which(duplicated(rownames(exprs.GBM))),]
exprs.GBM[1:5,1:5]

# list of genes
recur.mut.gene <- c("KRAS", "TP53", "STK11", "RBM10", "SPATA31C1", "KRTAP4-11",
                    "DCAF8L2", "AGAP6", "KEAP1", "SETD2", "ZNF679", "FSCB",
                    "BRAF", "ZNF770", "U2AF1", "SMARCA4", "HRNR", "EGFR")

nk.signature = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1")
set_A = c("NT5E")


# KEGG_hsa
load(system.file("data", "KEGG_hsa.rda", package = "gsean"))

pathways  = gmtPathways("/depot/pccr/data/Matosevic/tcga_analysis/Analysis/1_gsean/msigdb/NK_signature_custom.gmt")
pathways  = gmtPathways("/depot/pccr/data/Matosevic/tcga_analysis/Analysis/1_gsean/msigdb/c2.cp.kegg.v7.1.symbols.gmt")


# GSEA
set.seed(1)

result.GSEA <- gsean(KEGG_hsa, recur.mut.gene, exprs.GBM, threshold = 0.7)

result.GSEA.NK <- gsean(pathways, set_A, exprs.GBM, threshold = 0.7)
head(result.GSEA.NK)


fwrite(result.GSEA, file = "results.GSEA.txt", sep = "\t", sep2=c("", " ", ""))

result.GSEA %>%
  filter(str_detect(pathway, "NK"))

invisible(capture.output(p <- GSEA.barplot(result.GSEA, category = 'pathway',
                                           score = 'NES', pvalue = 'padj',
                                           sort = 'padj', top = 20)))



p <- GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES',
                  pvalue = 'padj', sort = 'padj', top = 20)
p + theme(plot.margin = margin(10, 10, 10, 75))

