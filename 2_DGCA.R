library(DGCA)
library(TCGAbiolinks)
library(EnsDb.Hsapiens.v86)
library(corrplot)
library(tidyverse)
library(RColorBrewer)

setwd("/depot/pccr/data/Matosevic/tcga_analysis/Analysis/2_DGCA")

query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM")

GDCdownload(query)
data <- GDCprepare(query,
                   save = TRUE, 
                   save.filename = "TCGA_GBM_HTSeq_FPKM.rda")

load(file = "TCGA_GBM_HTSeq_FPKM.rda")

samplesDown <- getResults(query,cols = c("cases"))
length(samplesDown)

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

length(dataSmTP)

dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")

length(dataSmNT)

dataPrep1 = data
dim(dataPrep1)


dataPrep <- TCGAanalyze_Preprocessing(object = dataPrep1, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - FPKM")

# dim(dataPrep)
# 
# dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                       geneInfo = geneInfoHT,
#                                       method = "gcContent")
# dim(dataNorm)
# 
# dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                       geneInfo = geneInfoHT,
#                                       method = "geneLength")
# dim(dataNorm)
# 
# dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
#                                   method = "quantile", 
#                                   qnt.cut =  0.25)
# 
# dim(dataFilt)
# 
# dataPrep_FPKM <- UseRaw_afterFilter(dataPrep, dataFilt)
# dim(dataPrep_FPKM)
# head(dataPrep_FPKM)

#tumor.FPKM  = dataPrep_FPKM[,dataSmTP]
#normal.FPKM = dataPrep_FPKM[,dataSmNT]


tumor.FPKM  = dataPrep[,dataSmTP]
normal.FPKM = dataPrep[,dataSmNT]

dim(tumor.FPKM)
dim(normal.FPKM)

saveRDS(tumor.FPKM, file  = "tumor.GBM.FPKM.rds")
saveRDS(normal.FPKM, file = "normal.GBM.FPKM.rds")

##############################################################################

tumor.FPKM = readRDS(file  = "tumor.GBM.FPKM.rds")
normal.FPKM = readRDS(file = "normal.GBM.FPKM.rds")
  
dim(tumor.FPKM)
dim(normal.FPKM)

head(normal.FPKM)

#Prepare design matrix
m1 = cbind(c(rep(1, 5)), c(rep(0, 5))  )
m2 = cbind(c(rep(0, 156)), c(rep(1, 156)))
m = rbind(m1,m2)


rownames(m) = c(colnames(normal.FPKM), colnames(tumor.FPKM))
colnames(m) = c("Normal", "Tumor")
head(m)
dim(m)

#prepare combined FPKM matrix
mat = cbind(normal.FPKM, tumor.FPKM)
mat = tumor.FPKM
dim(mat)
mat[1:5,1:5]

mat = as.data.frame(mat)
mat = rownames_to_column(mat, "Gene_ID")

ens2symbol <- AnnotationDbi::select(EnsDb.Hsapiens.v86,
                                    key=mat$Gene_ID, 
                                    columns="SYMBOL",
                                    keytype="GENEID")


ens2symbol <- as_tibble(ens2symbol)
names(ens2symbol) = c("Gene_ID", "SYMBOL")

ens2symbol %>%
  select(everything()) %>%
  summarise_all(funs(sum(is.na(.))))

x = inner_join(as.data.frame(mat), ens2symbol)
x = remove_rownames(x)
dim(x)


x = x %>%  
  drop_na(SYMBOL)

dim(x)

x =  dplyr::distinct(x, SYMBOL, .keep_all = TRUE)
x = column_to_rownames(x, "SYMBOL")
x = dplyr::select(x, -Gene_ID)
x[1:5, 1:4]


M = t(x)
M = as.data.frame(M)
M[1:5,1:5]
dim(M)
saveRDS(M, file  = "Tumor.Matrix.rds")


NKa = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "NT5E")
filename = "NKa.png"
data_sub = dplyr::select(M, NKa)
write.table(data_sub, file = "NKa.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()



NKb = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "HAVCR2")
filename = "NKb.png"
data_sub = dplyr::select(M, NKb)
write.table(data_sub, file = "NKb.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()


NKc = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "B4GALNT1")
filename = "NKc.png"
data_sub = dplyr::select(M, NKc)
write.table(data_sub, file = "NKc.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()


NKd = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "PVR")
filename = "NKd.png"
data_sub = dplyr::select(M, NKd)
write.table(data_sub, file = "NKd.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()


NKe = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "MICA", "MICB")
filename = "NKe.png"
data_sub = dplyr::select(M, NKe)
write.table(data_sub, file = "NKe.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()

NKf = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "B4GALNT1", "MICA", "MICB", "NT5E")
filename = "NKf.png"
data_sub = dplyr::select(M, NKf)
write.table(data_sub, file = "NKf.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()

NKg = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "IL15")
filename = "NKg.png"
data_sub = dplyr::select(M, NKg)
write.table(data_sub, file = "NKg.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()


NKh = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "CXCL10")
filename = "NKh.png"
data_sub = dplyr::select(M, NKh)
write.table(data_sub, file = "NKh.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()


NKi = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1", "CCL5")
filename = "NKi.png"
data_sub = dplyr::select(M, NKi)
write.table(data_sub, file = "NKi.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()


NKj = c("NCR1", "NCR3", "KLRB1", "CD160", "PRF1","NT5E", "HAVCR2", "B4GALNT1", "PVR", "MICA", "MICB", "IL15", "CXCL10", "CCL5")
filename = "NKj.png"
data_sub = dplyr::select(M, NKj)
write.table(data_sub, file = "NKj.txt", sep = "\t", quote = F)

data_sub_corr <- cor(data_sub)
png(filename, height = 12, width = 12, units = 'in', res = 300, pointsize = 12)
corrplot(data_sub_corr, method="color", tl.cex = 2, cl.cex = 1.5)
dev.off()

corrplot.mixed(data_sub_corr, tl.cex = 2, cl.cex = 1.)


corrplot(data_sub_corr, method="circle", tl.cex = 2, cl.cex = 1.5, 
         order = "hclust", addrect = 3,
         col = brewer.pal(n = 8, name = "RdBu"))




##############################################################################
library(ComplexHeatmap)
library(circlize)


#correlation Analysis
ddcor_res = ddcorAll(inputMat = mat, design = m,
                     compare = c("Tumor", "Normal"),
                     adjust = "none", nPerm = 0, nPairs = 100)

saveRDS(ddcor_res, file = "ddcor_res.rds")

head(ddcor_res)

plotVals(inputMat = mat, design = m,
         compare = c("Tumor", "Normal"),
         gene = "ENSG00000075624")

plotCors(inputMat = mat, design = m,
         compare = c("Tumor", "Normal"),
         geneA = "ENSG00000075624", geneB = "ENSG00000164796")
