# LASSO

# packages
library(tidyverse)
library(survival)
library(survminer)
library(randomForestSRC)
library(glmnet)

# data preparation ####


#load data
# gene expression
# downloaded from:
# https://xenabrowser.net/datapages/?cohort=TCGA%20Prostate%20Cancer%20(PRAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# gene expression RNAseq: IlluminaHiSeq (n=550) TCGA Hub
gene.expression.data <- read.delim("TCGA.PRAD.sampleMap_HiSeqV2.gz", row.names=1, check.names = F)

# selection of reduced number of genes for example
selected.genes=c( "ALDH5A1",
                 "CASP9",
                  "CCT3",
                  "CD24", 
                  "CDC25B",
                  "DNMT3B",
                 "DPPA",
                 "IL1RAP",
                  "IRF5",
                  "MEN1",
                 "MICAL3",
                 "RND3", 
                 "RPSA",
                 "KMT5C",
                 "TYMS", 
                 "ZMYND8")

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

df$PFI.time=df$PFI.time/30 #set time in months

# lasso ####

alpha <- 1
set.seed(777)

cv.fit <- cv.glmnet(as.matrix(df[selected.genes]), Surv(df$PFI.time, df$PFI), family="cox",  maxit = 1000, alpha=alpha)
coeffs.lasso <- coef(cv.fit, s = cv.fit$lambda.min)
save.names=coeffs.lasso@Dimnames[[1]]
coeffs.lasso <- as.numeric(coeffs.lasso)
names(coeffs.lasso)=save.names
coeffs.lasso

# evaluate model ####


mat=df[names(coeffs.lasso)]
# calculate score
for (i in seq_along(coeffs.lasso)) {
  col_name <- names(coeffs.lasso)[i]
  mat[[col_name]] <- mat[[col_name]] * coeffs.lasso[i]
}
mat$score=rowSums(mat)
df$score=mat$score

# assess score as a continuous variable
res.cox <- coxph(Surv(PFI.time, PFI) ~ score, data = df)
res.cox

# assess score as a dicotomized variable
res.cox <- coxph(Surv(PFI.time, PFI) ~ score>median(score, na.rm = T), data = df)
res.cox


## plot ####
cut=median(df$score, na.rm = T)
df$expression=ifelse(df$score>cut, "High", "Low") %>% factor(levels = c("Low", "High"))
table(df$expresion)
fit=survfit(Surv(PFI.time, PFI) ~ expression, data = df)

hr<- coxph(Surv(PFI.time, PFI) ~ expression, data = df) # calculate hazard ratio
summary(hr)

plot=ggsurvplot(fit,palette = "Dark2")
plot

plot=ggsurvplot(fit,
                ggtheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 25)),
                legend='none',
                title = expression(italic("TCGA-PRAD Risk Score")),
                legend.labs = c("low", "high"),
                conf.int = FALSE,
                censor.shape=124,
                censor.size=3,
                palette = c("#2DACA3","#D55AAA"),
                xlab='Months',
                ylab='Progression-Free Survival',
                font.y=c(15),
                break.time.by=12,
                font.tickslab=c(12),
                font.x=15,
                pval=T,
                tables.y.text = FALSE,
                risk.table = TRUE,
                cumevents = FALSE,
                tables.title= FALSE,
                tables.height = 0.15,
                tables.theme = theme_cleantable(),
                risk.table.title=element_blank()
                
)

plot
