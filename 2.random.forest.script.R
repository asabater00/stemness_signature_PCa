# RANDOM FOREST

# packages
library(tidyverse)
library(survival)
library(survminer)
library(randomForestSRC)

# data preparation ####

#load data
# gene expression
# downloaded from:
# https://xenabrowser.net/datapages/?cohort=TCGA%20Prostate%20Cancer%20(PRAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# gene expression RNAseq: IlluminaHiSeq (n=550) TCGA Hub
gene.expression.data <- read.delim("TCGA.PRAD.sampleMap_HiSeqV2.gz", row.names=1, check.names = F)

selected.genes=c("ABCB1", "ABCG2", "ACTB", "ADAM10", "ADAM15", "ADD1", "ADRM1", "ALDH1", "ALDH5A1", "ALDH7A1", "APMAP", "ARPC1B", "ATP2A3", "B4GALT4", "BCL2L1", "BMI", "BRCA1", "C1QBP", "CAP1", "CAPG", "CASP9", "CCNH", "CCT3", "CD133", "CD164", "CD24", "CD26", "CD44", "CD63", "CD9", "CDC25B", "CDON", "CETN2", "CHGA", "CIZ1", "KIT", "CNOT2", "CPT1B", "CREB3L1", "CTNNB1", "CTSB", "CXCR4", "DBNL", "DDX3X", "DGCR6", "DNMT3B", "DPPA", "EIF4A1", "FBLN1", "FBXO11", "FHL1", "FOXO4", "FXYD5", "GIPC1", "GNAI", "GNAS", "GPR125", "HDAC1", "HRH1", "HSPA5", "IKBKG", "IL1RAP", "ILF3", "IRF5", "ITGA6", "JAK2", "KDM4A", "KIFAP3", "KLF12", "KTR18", "KTR5", "LMNA", "LTBP1", "LTBR", "MAD1L1", "MAGED2", "MAX", "MB21D1", "MBNL2", "MDM2", "MEN1", "MICAL3", "MLX", "MYADM", "NANOG", "NARS", "NES", "NFATC3", "OCIAD1", "PFKFB3", "PLS3", "PMEPA1", "POU5F1B", "PRMT1", "PSEN1", "PTHLH", "PTPN18", "PTPRF", "RAPH1", "RASSF8", "RFC1", "RIT1", "RND3", "RNF145", "RPS3", "RPSA", "S6K1", "SAR1B", "SEMA3F", "SEMA4C", "SERPINE1", "SERPINE2", "SHISA2", "SIK1", "SIRT1", "SLC38A1", "SMARCA4", "SMO", "SOX2", "SOX9", "SPINT1", "STAT2", "STRAP", "STX1A", "SUPT5H", "KMT5C", "TAB3", "TACSTD2", "TAF4", "TAF5L", "TBCA", "TCF15", "TCF25", "TFCP2", "THOC5", "TNS3", "TP53INP2", "TRIM47", "TSC2", "TYMS", "UBE3C", "UBTD2", "VEGFA", "ZMYND8")

selected.genes=selected.genes[selected.genes%in%rownames(gene.expression.data)]

gex.data=t(gene.expression.data[selected.genes,])

gex.data=gex.data[complete.cases(gex.data),]

# pheno data

# TCGA pheno data selexted and downloaded from xena
# available at: https://xenabrowser.net/?bookmark=f304e43831b09a9db997926f3dc9ea27

data=read.delim("TCGA_from_xena.tsv", 
                stringsAsFactors=TRUE, 
                na.strings = c("", " ")) 

all(rownames(gex.data)%in%data$sample)

# clean data ####
data=data %>% 
  filter(sample_type=="Primary Tumor")

data$ISUP=case_when(data$gleason_score==6~1,
                    data$gleason_score==7 & data$primary_pattern==3~2,
                    data$gleason_score==7 & data$primary_pattern==4~3,
                    data$gleason_score==8~4,
                    data$gleason_score%in%c(9,10)~5,
)
data=data %>% column_to_rownames("sample")


df=merge(data, as.data.frame(gex.data), by=0) %>% column_to_rownames("Row.names")

df=df[,c("PFI.time","PFI",selected.genes)]

# random forest ####

set.seed(1)
tune=tune(Surv(PFI.time, PFI)~ ., data=df, verbose=1)
ns=as.numeric(tune$optimal[1])
print(ns)
mt=as.numeric(tune$optimal[2])
print(mt)

fit=rfsrc(Surv(PFI.time, PFI)~ ., data=df, mtry = mt, nodesize = ns, verbose=1, importance = T)


## see variable importance of genes ####
fit.vimp=fit$importance
order.vimp=names(sort(fit.vimp, decreasing = T))
res=list(importancia=fit.vimp, rank.importancia=data.frame(gen=order.vimp,ranking=c(1:length(order.vimp))))
res
