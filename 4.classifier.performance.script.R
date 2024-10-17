# packages
library(tidyverse)
library(gplots)
library(ggplot2)
library(pheatmap)
library(factoextra)
library(dplyr)
library(canvasXpress)
library(RColorBrewer)
library(pROC)

# download clinical and genomic data from cBioportal (https://www.cbioportal.org/study/summary?id=prad_msk_mdanderson_2023). 
# The same script was used (with adaptations) for the Beltran et al. dataset analysis (https://www.cbioportal.org/study/summary?id=prad_fhcrc).

# set working directory where files from cBioPortal are
setwd("")

# load data ####
exp<-read.delim("prad_msk_mdanderson_2023/prad_msk_mdanderson_2023/data_mrna_seq_v2_rsem.txt")
pheno<-read.delim("prad_msk_mdanderson_2023/prad_msk_mdanderson_2023/data_clinical_sample.txt")

# clean data ####
pheno<-pheno[-c(1:4),]
pheno<-pheno[-grep("T200",pheno$X.Sample.Identifier),]
colnames(pheno)[1]<-"ID"

exp<-exp[,c(1:46)]
exp<-exp[-which(duplicated(exp$Hugo_Symbol)),]
rownames(exp)<-exp$Hugo_Symbol

exp<-exp[,-c(1,2)]
exp<-as.data.frame(t(exp))
exp<-exp %>% mutate (ID=rownames(exp)) %>% relocate(ID)

exp$ID <- gsub("\\.", "-", exp$ID)

# merge expression and phenotipic data
data<-merge(pheno, exp, by="ID")

# select genes of interest
genes<-c("KMT5C", "DNMT3B", "IRF5", "DPP4", "MEN1", "TYMS", "CDC25B")

# filter dataframe
data<-data[,c(1:12,which(colnames(data)%in%genes))]

# UNSUPERVISED CLUSTERING ####

# calculate zscore
mat_zscore <- scale(data [,c(13:19)])

range <- max(abs(mat_zscore))
range

# choose columns for annotations
annot.heatmap <- data.frame((data[,c(6,7)]))
rownames(annot.heatmap) <- rownames(mat_zscore)<-data[,1]

## PLOT ####

sampletype <- c("Primary"="#0096c7", "Metastasis"="#fdd128")
annot.heatmap$Sample.Type<-factor(annot.heatmap$Sample.Type, levels=c("Primary", "Metastasis"))
morphology <- c("Adenocarcinoma"="#db0f27", "Sarcomatoid Carcinoma"="burlywood3", "Neuroendocrine Carcinoma"= "#741791", "Poorly Differentiated Carcinoma"="#4ea697")
annot.heatmap$Morphology<-factor(annot.heatmap$Morphology, levels=c("Adenocarcinoma", "Sarcomatoid Carcinoma", "Neuroendocrine Carcinoma", "Poorly Differentiated Carcinoma"))

annotation_colors <- list(Sample.Type=sampletype, Morphology=morphology)

heat <- pheatmap (t(mat_zscore), 
                  color =colorRampPalette(c("dodgerblue3","white", "firebrick3"))(100)[1:100],
                  border_color = "gray", 
                  fontsize = 7,
                  show_rownames = TRUE, 
                  cluster_cols = TRUE, 
                  cluster_rows = TRUE,
                  breaks = seq(-range, range, length.out = 100),
                  annotation_col = annot.heatmap,
                  annotation_colors = annotation_colors, 
                  annotation_legend = T)

heat

# PRINCIPAL COMPONENT ANALYSIS ####

# prepare table
table = data
rownames(table) <- table[,1]
table2 <- table

# choose genes
my_data <- table2[,13:19]

res.pca <- prcomp(my_data, scale = TRUE)
res.pca

table2$Morphology<-factor(table2$Morphology, levels=c("Adenocarcinoma", "Sarcomatoid Carcinoma", "Neuroendocrine Carcinoma", "Poorly Differentiated Carcinoma"))

## PLOT ####

p<- fviz_pca_biplot(res.pca, habillage = table2$Morphology,
                    geom = "point", pointsize=2,
                    label = "var",
                    addEllipses = FALSE, 
                    palette = c('#db0f27', 'burlywood3', '#741791', "#4ea697"), 
                    geom.var = '',
                    text = element_text(size=54),
                    
)+
  scale_shape_manual(values=c(19,19,19,19))

p

Contrib <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 15, fill="gray", color = "black")

Contrib

# RECEIVER OPERATING CHARACTERISTICS (ROC) CURVE ####

data$Morphology %>% table()

data$group=ifelse(data$Morphology=="Neuroendocrine Carcinoma",1,0)
data$group=as.factor(data$group)

# calculate risk score
data$Risk.score<- 0.284 * data$KMT5C - 0.060 * data$DPP4 + 0.218 * data$TYMS + 0.048 * data$CDC25B + 0.090 * data$IRF5 + 0.272 * data$MEN1 + 0.083 * data$DNMT3B

# model
roc.model=roc(data$group ~ data$Risk.score)
auc_value <- auc(roc.model)
auc_ci <- ci.auc(roc.model)

## PLOT ####
roc_data <- data.frame(
  specificity = rev(roc.model$specificities),  # Reverse the order for ggplot
  sensitivity = rev(roc.model$sensitivities)
)

p<-ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#741791", linewidth = 1) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "darkgray") + 
  labs(
    title = paste0("7-gene score and NEPC PDXs"),
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme(plot.title=element_text(hjust=0.75))+
  scale_x_reverse()+
  theme_minimal() + 
  coord_equal()+
  annotate("text", x = 0.1, y = 0.1, label = paste("AUC =", round(auc_value,2), "(95%CI = 0.84-1)"), size = 7, hjust = 0.9, vjust = 0) + 
  theme(
    plot.title = element_text(hjust = 0.5, size=20),
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20),
    axis.text.x=element_text(size = 14),
    axis.text.y=element_text(size = 14)
  ) 

p
