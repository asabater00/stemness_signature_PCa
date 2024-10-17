# SURVIVAL ANALYSIS

# packages
library(tidyverse)
library(survival)
library(survminer)
library(pheatmap)


# load data ####

# gene expression
# downloaded from:
# https://xenabrowser.net/datapages/?cohort=TCGA%20Prostate%20Cancer%20(PRAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# gene expression RNAseq: IlluminaHiSeq (n=550) TCGA Hub

gene.expression.data <- read.delim("TCGA.PRAD.sampleMap_HiSeqV2.gz", row.names=1, check.names = F)

# selected genes for this example
selected.genes=c("DBNL", "UBTD2", "MBNL2")
gex.data=t(gene.expression.data[selected.genes,])
gex.data=gex.data[complete.cases(gex.data),]

# pheno data

# TCGA pheno data selected and downloaded from XENA
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

# survival analysis ####
# using DBNL for this example
gene.of.interest="DBNL"


# cut-off point is calculated using CutoffFinder (DOI: 10.1371/journal.pone.0051862)
cut.point=11.1715
table(df$gene.expression)

## univariable survival ####
df=df %>% 
  mutate(gene.expression=factor(ifelse(get(gene.of.interest)>cut.point,1,0)))

# cox analysis
coxph(Surv(time = PFI.time, event = PFI)~gene.expression, df)

### PLOT KM ####
fit=survfit(Surv(time = PFI.time, event = PFI)~gene.expression, df)
ggsurvplot(fit)


## multivariable ####

coxph(Surv(time = PFI.time, event = PFI)~gene.expression + ISUP + targeted_molecular_therapy + radiation_therapy , df)
