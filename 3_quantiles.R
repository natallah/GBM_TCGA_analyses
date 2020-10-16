##############################################################################
#       Quantile normalization
##############################################################################
library(tidyverse)
library(TCGAbiolinks)
library(EnsDb.Hsapiens.v86)
library(fgsea)
library(VennDiagram)
library(gplots)
library(ggplot2)

setwd("/depot/pccr/data/Matosevic/tcga_analysis/Analysis/3_quantile")
##############################################################################
#    				Functions    
##############################################################################
# We have followiung three functions"

# 1 - classify_patients
# This function takes the gene as input and classify the patients into high and low groups. 
# The results is a list with 4 variables as below
# counts - counts of patients by quantile
# summary - (min, max, mean, median) by quantile
# high - patient IDs with high expression of given gene
# low  - patient IDss with low expression of given gene

# 2 - calculate_DEG
# This function takes treatment and control groups (in order) as input i.e. $high and $low from function classify_patients. 
# The results is a dataframe of DE results containing columns (Gene_ID, logFC, logCPM, PValue, FDR, SYMBOL)


# 3 - perform_fgsea
# This function takes DEG table (output from function calculate_DEG) and pathways (gmt file from MSigDB) as input and performs FGSEA analsysi. 
# The results is FGSEA enrichment results performed with given DEG table and pathways of interest.

classify_patients <- function(gene) {
  
  gene_data =  dplyr::select(M, gene)
  gene_data =  rownames_to_column(gene_data, "patient_ID")
  names(gene_data) = c("patient_ID", "gene")
  head(gene_data)
  
  #print(mean(gene_data$gene))
  #print(median(gene_data$gene))
  
  data_10Q = gene_data %>% 
    mutate(quintile = ntile(gene, 10))
  
  counts = data_10Q %>% 
    group_by(quintile)  %>% 
    summarize(n())
  
  summary = data_10Q %>% 
    group_by(quintile)  %>% 
    summarize(size_min = min(gene), size_mean = mean(gene), sie_median = median(gene), size_max = max(gene))
  
  high = data_10Q %>%
    dplyr::filter( quintile > 5) %>%
    dplyr::select("patient_ID")
  
  low  = data_10Q %>%
    dplyr::filter( quintile <= 5) %>%
    dplyr::select("patient_ID")
  
  high = as.vector(high$patient_ID)
  low  = as.vector(low$patient_ID)
  
  result <- list(counts=counts,summary=summary, high = high, low = low)
  return(result)
  
}

calculate_DEG <- function(treatment, control) {
  
  # extract matrix of raw counts
  treatment  = dataPrep_raw[,treatment]
  control   = dataPrep_raw[,control]
  
  print(lapply(list(treatment, control), dim))
  
  #calculate deg as high_vs_low
  
  deg = TCGAanalyze_DEA(
    mat1 = control,
    mat2 = treatment,
    Cond1type = "control",
    Cond2type = "treatment",
    fdr.cut = 1
  )
  
  deg = rownames_to_column(deg, "Gene_ID")
  
  #Link Gene Symbol Information
  ens2symbol <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                      key=deg$Gene_ID, 
                                      columns="SYMBOL",
                                      keytype="GENEID")
  
  ens2symbol <- as_tibble(ens2symbol)
  names(ens2symbol) = c("Gene_ID", "SYMBOL")
  deg <- inner_join(deg, ens2symbol)
  
  return(deg)
  
}

perform_fgsea <- function(deg_table, pathways) {
  
  #calculate Ranks
  deg_table$RANK = sign(deg_table$logFC) * -log10(deg_table$PValue)
  
  #calculate GSEA Ranks
  GSEA <- deg_table %>% 
    dplyr::select(SYMBOL, RANK) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL)  %>% 
    summarize(RANK=max(RANK))
  
  GSEA = arrange(GSEA, desc(RANK))
  GSEA.ranks <- deframe(GSEA)
  head(GSEA.ranks)
  
  set.seed(1)
  GSEA.results <- fgsea(pathways=pathways, stats=GSEA.ranks, nperm=1000)
  GSEA.results$direction = if_else(GSEA.results$NES >= 0, "Positive", "Negative")
  GSEA.results = arrange(GSEA.results, pval)
  
  return(GSEA.results)
  
}

plot_fgsea <- function(GSEA.results) {

	fgsea.sig.pos = fgsea.results %>%
  	filter(pval <= 0.05, direction == "Positive") %>%
  	top_n(-10, pval)

	fgsea.sig.neg = fgsea.results %>%
  	filter(pval <= 0.05, direction == "Negative") %>%
  	top_n(-10, pval)

	fgsea.sig = rbind(fgsea.sig.pos, fgsea.sig.neg)
	
	p = ggplot(fgsea.sig, aes(x = reorder(pathway, NES), y = NES)) +
  	geom_col(aes(fill = direction)) +
  	coord_flip( ylim = c(-2,2) ) +
  	labs(x="Pathway", y="Normalized Enrichment Score",
        title="Significant Gene Sets", fill = "NES") + 
  	theme_minimal()

	p + scale_fill_manual(values=c("#999999", "#E69F00"))

	return(p)

}


##############################################################################
#    				Read_Data    
##############################################################################


# M is the tumor.FPKM matrix for GBM patients
# This matrix is used for clsiifying the patients into high and low groups.
# FPKM matrix is already normalized and better than using the raw counts.
M = readRDS(file  = "Tumor.Matrix.rds")
M[1:5,1:5]

# dataPrep_raw is the matrix of raw counts from GBM data.
# The raw counts are used for DEG analysis which performs internal normalization.
dataPrep_raw = readRDS(file  = "dataPrep_raw.rds")

# pre-computed names for the 5 GBM normal samples
# This is used in case we need to compare High/low against the normals
GBM_normal = c("TCGA-06-0675-11A-32R-A36H-07", "TCGA-06-0678-11A-32R-A36H-07", 
               "TCGA-06-0681-11A-41R-A36H-07", "TCGA-06-0680-11A-32R-A36H-07", 
               "TCGA-06-AABW-11A-31R-A36H-07")

##############################################################################
#    				Analysis    
##############################################################################

gene = "MICB"
high_out = paste0(gene,"_high.venn")
low_out = paste0(gene,"_low.venn")

p = classify_patients(gene)
write.table(p$high, file = high_out, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(p$low,  file = low_out, sep = "\t", quote = F, row.names = F, col.names = F)
 

p = classify_patients("NT5E")
q = classify_patients("HAVCR2")

high_vs_low = calculate_DEG(p$high, p$low)

pathways  = gmtPathways("/depot/pccr/data/Matosevic/tcga_analysis/Analysis/1_gsean/msigdb/NK_signature_custom.gmt")
pathways  = gmtPathways("/depot/pccr/data/Matosevic/tcga_analysis/Analysis/1_gsean/msigdb/c2.cp.kegg.v7.1.symbols.gmt")
gsea.results = perform_fgsea(high_vs_low, pathways)
gsea.results[1:10,1:5]


#Venn
mylist = list(p$low, q$low)
mylist = list(p$high, q$high)

ItemsList  <- venn(mylist, show.plot=FALSE)
x = attr(ItemsList,"intersections")

low.intersec = x[["A:B"]]
high.intersec = x[["A:B"]]

intesec_list = list(low.intersec, high.intersec)
  
venn.diagram(mylist, 
             fill = c("#FEAD72", "#FED976"),
             category.names = c("NT5E_high", "HAVCR2_high"),
             output = TRUE ,
             imagetype="png" ,
             height = 2000 , 
             width = 2000 , 
             resolution = 300,
             filename = "Venn.png",
             main = "Venn",
             main.cex = 1.5, main.just = c(1,-2),
             alpha = c(0.5),
             scaled = FALSE, euler.d = FALSE, # This is needed if sets are inclusive
             cat.pos  = c(0, 180),
             #cat.just     = list(c(0, -0.7), c(0, 0)),
             cex = 1.5, cat.cex = 1.5,
             print.mode = c("raw", "percent")
)

