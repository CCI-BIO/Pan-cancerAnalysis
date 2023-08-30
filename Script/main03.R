# Local
# Clear the environment
rm(list = ls())

# Load necessary libraries if any
library(stringr)
library(dplyr)
require(data.table)
library(readxl)
library(miRLAB)
library(miRBaseConverter)
library(ggplot2)
library(varhandle)
library(scales)
library(reshape)
library(plyr)
library(RColorBrewer)
library(tidyverse)
library(xtable)
library(TCGAbiolinks)
library("optparse")
library(grid)
library(futile.logger)
library(VennDiagram)
library(ezcox)
library(glmnet)
library(survival)
library(survminer)
library(timeROC)

#---------------------------------------
# Set environment variables if any
# Please remember to create necessary folders
# LGG
# dir.create(file.path("C:/Users/vpham/Documents/002NetworkAnalysis/Data", "LGG"), showWarnings = FALSE)
# rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# outDir="C:/Users/vpham/Documents/002NetworkAnalysis/Data/LGG"
rootDir="/srv/scratch/z3538133/002NetworkAnalysis" # And put the input files in "rootDir/Data"
cancertype="LGG"
# cancertype="GBM"
outDir=paste(rootDir, "/Data/", cancertype, sep = "")
# controlDir="C:/Users/vpham/Documents/002NetworkAnalysis/control"
# primarySite=""
#---------------------------------------

# Set working directory
setwd(rootDir)

# Include the script of functions
source(paste(rootDir, "/Script/ProposedMethod_Functions.R", sep=""))

#================================================================
# (5) Constructing and validating risk-score system
#================================================================

######constructing risk-score system, LASSO model#######

load(paste(outDir, '/coxdata.Rdata', sep = ""))
# > coxdata[1:5,1:6]
# OS OS.time    ALB  ATP7A   CASP3    CDK1
# TCGA.CS.4938.01  0    3574 5.4263 8.4838 10.3608  6.3923
# TCGA.CS.4941.01  1     234 5.7808 9.1396 10.8556  9.0084
# TCGA.CS.4942.01  1    1335 4.0875 9.3772 10.6027  8.9425
# TCGA.CS.4943.01  1    1106 3.8074 9.5018 10.4798 12.3250
# TCGA.CS.4944.01  0    1828 2.3219 7.8765  9.8281  3.9069
# > nrow(coxdata)
# [1] 348
# > ncol(coxdata)
# [1] 15
# > colnames(coxdata)
# [1] "OS"       "OS.time"  "ALB"      "ATP7A"    "CASP3"    "CDK1"     "CP"       "CYP1A1"   "F5"       "MAP1LC3A" "MT-CO1"   "SNCA"     "SP1"     
# [14] "TP53"     "XIAP"    

#=====================
# Only get genes from Univariate Cox analysis, from file "LGG_prognostic_genes.csv"
uni_gene <- read.csv(file = paste(outDir, "/", cancertype, "_prognostic_genes.csv", sep =""))
newList <- uni_gene$x
coxdata <- coxdata[,c(1,2,which(colnames(coxdata) %in% newList))]
coxdata[1:5,1:4]
nrow(coxdata)
ncol(coxdata)
# > coxdata[1:5,1:4]
# OS OS.time    ALB  ATP7A
# TCGA.CS.4938.01  0    3574 5.4263 8.4838
# TCGA.CS.4941.01  1     234 5.7808 9.1396
# TCGA.CS.4942.01  1    1335 4.0875 9.3772
# TCGA.CS.4943.01  1    1106 3.8074 9.5018
# TCGA.CS.4944.01  0    1828 2.3219 7.8765
# > nrow(coxdata)
# [1] 348
# > ncol(coxdata)
# [1] 12

print(paste(uni_gene$x, collapse=", "))
# > print(paste(uni_gene$x, collapse=", "))
# [1] "ALB, CASP3, CDK1, ATP7A, MT-CO1, CP, CYP1A1, F5, SP1, TP53"
#================

ncol <- ncol(coxdata)
irg_expr <- coxdata[ ,c(3:ncol)]

x <- as.matrix(irg_expr)
y <- data.matrix(Surv(coxdata$OS.time,coxdata$OS))
x[1:5, 1:4]
nrow(x)
head(y)
nrow(y)
# > x[1:5, 1:4]
# ALB  ATP7A   CASP3    CDK1
# TCGA.CS.4938.01 5.4263 8.4838 10.3608  6.3923
# TCGA.CS.4941.01 5.7808 9.1396 10.8556  9.0084
# TCGA.CS.4942.01 4.0875 9.3772 10.6027  8.9425
# TCGA.CS.4943.01 3.8074 9.5018 10.4798 12.3250
# TCGA.CS.4944.01 2.3219 7.8765  9.8281  3.9069
# > nrow(x)
# [1] 348
# > head(y)
# time status
# [1,] 3574      0
# [2,]  234      1
# [3,] 1335      1
# [4,] 1106      1
# [5,] 1828      0
# [6,]   NA      0
# > nrow(y)
# [1] 348

# remove time with NA or == 0
idx <- which(is.na(y[,1]) | y[,1] <= 0)
if(length(idx) > 0) {
  x <- x[-idx,]
  y <- y[-idx,]
}
idx
x[1:5, 1:4]
nrow(x)
head(y)
nrow(y)
# > idx
# [1]   6 249 250
# > x[1:5, 1:4]
# ALB  ATP7A   CASP3    CDK1
# TCGA.CS.4938.01 5.4263 8.4838 10.3608  6.3923
# TCGA.CS.4941.01 5.7808 9.1396 10.8556  9.0084
# TCGA.CS.4942.01 4.0875 9.3772 10.6027  8.9425
# TCGA.CS.4943.01 3.8074 9.5018 10.4798 12.3250
# TCGA.CS.4944.01 2.3219 7.8765  9.8281  3.9069
# > nrow(x)
# [1] 345
# > head(y)
# time status
# [1,] 3574      0
# [2,]  234      1
# [3,] 1335      1
# [4,] 1106      1
# [5,] 1828      0
# [6,] 1222      0
# > nrow(y)
# [1] 345

# Coefficient profiles in the LASSO regression model.
pdf(file = paste(outDir, "/", cancertype, '_coefPro.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
fit0 <- glmnet(x, y, family = "cox", alpha = 1, nlambda = 1000)
plot(fit0)
dev.off()
pdf(file = paste(outDir, "/", cancertype, '_coefProLambda.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
plot(fit0, xvar="lambda", label=FALSE)
dev.off()

#==========================================
# Test data for LGG
exp_file <- paste(rootDir, "/Data/CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt", sep = "")
clinical_file <- paste(rootDir, "/Data/CGGA.mRNAseq_693_clinical.20200506.txt", sep = "")

new_ds_x <- fread(exp_file, fill = TRUE)
new_ds_y <- read.table(clinical_file, sep = "\t", fill = TRUE, header = T)
new_ds_x[1:5,1:6]
new_ds_y[1:5,1:6]
# > new_ds_x[1:5,1:6]
# gene_name CGGA_1002 CGGA_1003 CGGA_1010 CGGA_1012 CGGA_1014
# 1:      A1BG        21         8        13        26        18
# 2:  A1BG-AS1        36         3        19        51        46
# 3:      A1CF         3         2         3         4         1
# 4:       A2M      5968       800      1227     19094      4042
# 5:   A2M-AS1         7         4        18         6         9
# > new_ds_y[1:5,1:6]
# CGGA_ID  PRS_type Histology   Grade Gender Age
# 1 CGGA_1002   Primary        AA WHO III Female  43
# 2 CGGA_1003   Primary         O  WHO II Female  47
# 3 CGGA_1010   Primary         A  WHO II   Male  45
# 4 CGGA_1012 Recurrent        rO  WHO II   Male  45
# 5 CGGA_1014   Primary         A  WHO II   Male  42

new_ds_x <- as.data.frame(new_ds_x)
new_ds_y <- as.data.frame(new_ds_y)
rownames(new_ds_x) <- new_ds_x$gene_name
new_ds_x <- new_ds_x[,-1]
new_ds_x[1:5,1:6]
# > new_ds_x[1:5,1:6]
# CGGA_1002 CGGA_1003 CGGA_1010 CGGA_1012 CGGA_1014 CGGA_1017
# A1BG            21         8        13        26        18        26
# A1BG-AS1        36         3        19        51        46        48
# A1CF             3         2         3         4         1         4
# A2M           5968       800      1227     19094      4042      8167
# A2M-AS1          7         4        18         6         9        46

new_ds_x = log2(new_ds_x + 1)
new_ds_x[1:5,1:6]
# > new_ds_x[1:5,1:6]
# CGGA_1002 CGGA_1003 CGGA_1010 CGGA_1012 CGGA_1014 CGGA_1017
# A1BG      4.459432  3.169925  3.807355  4.754888  4.247928  4.754888
# A1BG-AS1  5.209453  2.000000  4.321928  5.700440  5.554589  5.614710
# A1CF      2.000000  1.584963  2.000000  2.321928  1.000000  2.321928
# A2M      12.543274  9.645658 10.262095 14.220907 11.981210 12.995767
# A2M-AS1   3.000000  2.321928  4.247928  2.807355  3.321928  5.554589

new_ds_x <- new_ds_x[which(rownames(new_ds_x) %in% colnames(x)),]
new_ds_x <- t(new_ds_x)
nGenes <- ncol(new_ds_x)
new_ds_x[1:5,1:6]
nrow(new_ds_x)
ncol(new_ds_x)
# > new_ds_x[1:5,1:6]
# ALB    ATP7A     CASP3     CDK1       CP   CYP1A1
# CGGA_1002 6.247928 0.000000  9.407268 8.867279 5.491853 3.584963
# CGGA_1003 6.266787 1.584963  8.299208 7.930737 3.584963 4.169925
# CGGA_1010 1.584963 7.451211  5.754888 6.781360 4.807355 0.000000
# CGGA_1012 4.321928 4.754888 10.227616 9.715962 6.459432 1.000000
# CGGA_1014 7.000000 3.169925 10.493855 8.826548 6.918863 3.700440
# > nrow(new_ds_x)
# [1] 693
# > ncol(new_ds_x)
# [1] 10

new_ds_y <- new_ds_y[, c("CGGA_ID","OS", "Censor..alive.0..dead.1.", "Histology")]
head(new_ds_y)
nrow(new_ds_y)
# > head(new_ds_y)
# CGGA_ID   OS Censor..alive.0..dead.1. Histology
# 1 CGGA_1002  305                        1        AA
# 2 CGGA_1003 3817                        0         O
# 3 CGGA_1010  246                        1         A
# 4 CGGA_1012 3679                        1        rO
# 5 CGGA_1014  263                        1         A
# 6 CGGA_1017  768                        1       GBM
# > nrow(new_ds_y)
# [1] 693

# convert rownames in counts dataframe into a column
new_ds_x <- as.data.frame(new_ds_x)
new_ds_x$sample <- rownames(new_ds_x)
new_ds_x <- new_ds_x %>% dplyr::select(sample, everything())
new_ds_x[1:5,1:6]
# > new_ds_x[1:5,1:6]
# sample      ALB    ATP7A     CASP3     CDK1       CP
# CGGA_1002 CGGA_1002 6.247928 0.000000  9.407268 8.867279 5.491853
# CGGA_1003 CGGA_1003 6.266787 1.584963  8.299208 7.930737 3.584963
# CGGA_1010 CGGA_1010 1.584963 7.451211  5.754888 6.781360 4.807355
# CGGA_1012 CGGA_1012 4.321928 4.754888 10.227616 9.715962 6.459432
# CGGA_1014 CGGA_1014 7.000000 3.169925 10.493855 8.826548 6.918863

data <- merge(new_ds_x, new_ds_y, by.x = "sample", by.y = "CGGA_ID")
# > colnames(data)
# [1] "sample"                   "ALB"
# [3] "ATP7A"                    "CASP3"
# [5] "CDK1"                     "CP"
# [7] "CYP1A1"                   "F5"
# [9] "MT-CO1"                   "SP1"
# [11] "TP53"                     "OS"
# [13] "Censor..alive.0..dead.1." "Histology"
data <- data %>% dplyr::select((nGenes+2):(nGenes+3), everything())
data <- data %>% dplyr::select(sample, everything())
data[1:5,1:6]
nrow(data)
ncol(data)
# > data[1:5,1:6]
# sample   OS Censor..alive.0..dead.1.      ALB    ATP7A     CASP3
# 1 CGGA_1002  305                        1 6.247928 0.000000  9.407268
# 2 CGGA_1003 3817                        0 6.266787 1.584963  8.299208
# 3 CGGA_1010  246                        1 1.584963 7.451211  5.754888
# 4 CGGA_1012 3679                        1 4.321928 4.754888 10.227616
# 5 CGGA_1014  263                        1 7.000000 3.169925 10.493855
# > nrow(data)
# [1] 693
# > ncol(data)
# [1] 14

# Update new_ds_x & new_ds_y
table(data$Histology)
# > table(data$Histology)
# 
# A   AA   AO  AOA  GBM    O   OA   rA  rAA  rAO rAOA rGBM   rO  rOA
# 85   82   46   16  140   45    8   34   70   36    5  109   15    1
if(cancertype == "LGG") {
  selectedId <- which(!(data$Histology %in% c(NA, "GBM", "rGBM")))
} else if(cancertype == "GBM") {
  selectedId <- which(data$Histology %in% c("GBM", "rGBM"))
}
data <- data[selectedId,]
data <- data[,-ncol(data)] # remove Histology column
data[1:5,1:6]
nrow(data)
ncol(data)
# > data[1:5,1:6]
# sample   OS Censor..alive.0..dead.1.      ALB    ATP7A     CASP3
# 1 CGGA_1002  305                        1 6.247928 0.000000  9.407268
# 2 CGGA_1003 3817                        0 6.266787 1.584963  8.299208
# 3 CGGA_1010  246                        1 1.584963 7.451211  5.754888
# 4 CGGA_1012 3679                        1 4.321928 4.754888 10.227616
# 5 CGGA_1014  263                        1 7.000000 3.169925 10.493855
# > nrow(data)
# [1] 443
# > ncol(data)
# [1] 13

# remove time with NA or <= 0
idx <- which(is.na(data[,2]) | data[,2] <= 0)
if(length(idx) > 0) {
  data <- data[-idx,]
}
idx
data[1:5, 1:4]
nrow(data)
# > idx
# [1]  18  21  24 156 207 215 218 258 264 283 298 319 324 325 338 340 344 356 358
# [20] 387 392 436 437
# > data[1:5, 1:4]
# sample   OS Censor..alive.0..dead.1.      ALB
# 1 CGGA_1002  305                        1 6.247928
# 2 CGGA_1003 3817                        0 6.266787
# 3 CGGA_1010  246                        1 1.584963
# 4 CGGA_1012 3679                        1 4.321928
# 5 CGGA_1014  263                        1 7.000000
# > nrow(data)
# [1] 420

#convert first column to rowname
rownames(data) <- data[,1]
data <- data[,-1]
colnames(data)[1:2] <- c("time", "status")
data[1:5,1:6]
# > data[1:5,1:6]
# time status      ALB    ATP7A     CASP3     CDK1
# CGGA_1002  305      1 6.247928 0.000000  9.407268 8.867279
# CGGA_1003 3817      0 6.266787 1.584963  8.299208 7.930737
# CGGA_1010  246      1 1.584963 7.451211  5.754888 6.781360
# CGGA_1012 3679      1 4.321928 4.754888 10.227616 9.715962
# CGGA_1014  263      1 7.000000 3.169925 10.493855 8.826548

new_ds_x <- as.matrix(data[,3:ncol(data)])
new_ds_y <- data.matrix(Surv(data[,1],data[,2]))
new_ds_x[1:5,1:6]
head(new_ds_y)
# > new_ds_x[1:5,1:6]
# ALB    ATP7A     CASP3     CDK1       CP   CYP1A1
# CGGA_1002 6.247928 0.000000  9.407268 8.867279 5.491853 3.584963
# CGGA_1003 6.266787 1.584963  8.299208 7.930737 3.584963 4.169925
# CGGA_1010 1.584963 7.451211  5.754888 6.781360 4.807355 0.000000
# CGGA_1012 4.321928 4.754888 10.227616 9.715962 6.459432 1.000000
# CGGA_1014 7.000000 3.169925 10.493855 8.826548 6.918863 3.700440
# > head(new_ds_y)
# time status
# [1,]  305      1
# [2,] 3817      0
# [3,]  246      1
# [4,] 3679      1
# [5,]  263      1
# [6,] 2527      1
#==========================================

#==========================================
# Update training dataset
if(length(colnames(x)) > length(colnames(new_ds_x))) {
  print("Not enough genes for testing")
  print("Missing genes")
  print(colnames(x)[which(!(colnames(x) %in% colnames(new_ds_x)))])
} else {
  print("Enough genes for testing")
}
x <- x[,which(colnames(x) %in% colnames(new_ds_x))]
x[1:5,1:6]
nrow(x)
ncol(x)
# > x[1:5,1:6]
# ALB  ATP7A   CASP3    CDK1      CP CYP1A1
# TCGA.CS.4938.01 5.4263 8.4838 10.3608  6.3923 11.5159 3.3219
# TCGA.CS.4941.01 5.7808 9.1396 10.8556  9.0084 10.2872 2.0000
# TCGA.CS.4942.01 4.0875 9.3772 10.6027  8.9425 10.1364 6.4757
# TCGA.CS.4943.01 3.8074 9.5018 10.4798 12.3250  8.8437 1.0000
# TCGA.CS.4944.01 2.3219 7.8765  9.8281  3.9069 13.0533 3.8074
# > nrow(x)
# [1] 345
# > ncol(x)
# [1] 10
#==========================================

set.seed(1)

# Cross-validation for tuning parameter screening in the LASSO regression model.
pdf(file = paste(outDir, "/", cancertype, '_crossVal.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
cv.fit <- cv.glmnet(x, y,
                    family="cox",
                    maxit = 1000000,
                    alpha=1)
print(cv.fit)
# > print(cv.fit)
# 
# Call:  cv.glmnet(x = x, y = y, family = "cox", maxit = 1e+06, alpha = 1) 
# 
# Measure: Partial Likelihood Deviance 
# 
# Lambda Index Measure     SE Nonzero
# min 0.02459    19   10.42 0.5443       6
# 1se 0.13121     1   10.83 0.4943       0
plot(cv.fit)
dev.off()

# LASSO_gene table
pdf(file = paste(outDir, "/", cancertype, '_lasso_fit.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
fit <- glmnet(x, y, alpha = 1, family='cox',lambda=cv.fit$lambda.min)
plot(fit)
coef(fit)
dev.off()

Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Coefficients

Active.Index
Active.Coefficients

lasso_gene <- row.names(Coefficients)[Active.Index]
lasso_min <- data.frame(Active.Index,Active.Coefficients,lasso_gene)
lasso_min
# > lasso_min
# Active.Index Active.Coefficients lasso_gene
# 1            1         -0.12157533        ALB
# 2            3          0.25051154      CASP3
# 3            4          0.17369917       CDK1
# 4            5          0.15084949         CP
# 5            6         -0.03929607     CYP1A1
# 6            8         -0.03095952     MT-CO1

save(lasso_min,file = paste(outDir, "/", cancertype, "_lasso_min.Rdata", sep = ""))
save(cv.fit,fit,lasso_gene,file = paste(outDir, "/", cancertype, "_Lasso_model_min.Rdata", sep = ""))

## validation of LASSO model ###### 
# NOTE: Should be run on a new dataset. new x and y should be assigned for new dataset.
# to validate with CGGA dataset with new gene expression which is the x variable 
lasso.prob <- predict(cv.fit, 
                      newx = as.matrix(new_ds_x),
                      s = c(cv.fit$lambda.min,cv.fit$lambda.1se))
re <- cbind(new_ds_y, lasso.prob)
head(re)
# > head(re)
# time status        1 2
# CGGA_1002  305      1 3.205046 0
# CGGA_1003 3817      0 2.401783 0
# CGGA_1010  246      1 2.757706 0
# CGGA_1012 3679      1 4.072335 0
# CGGA_1014  263      1 3.598046 0
# CGGA_1018 2527      1 2.845748 0

##### timeROC curve #####

data$Risk_score <- as.numeric(lasso.prob[ ,1])

with(data,
     ROC <<- timeROC(T = time,
                     delta = status,
                     marker = Risk_score,
                     cause = 1,
                     weighting = "marginal",
                     times = c(365, 1095, 1825),
                     ROC = TRUE,
                     iid = TRUE)
)

plot(ROC,time = 1825, col = "blue",add = FALSE)
plot(ROC,time = 1095,col = "blue",add = FALSE)
plot(ROC,time = 365,col = "blue",add = FALSE)

ROC$AUC
confint(ROC)

pdf(file = paste(outDir,"/", cancertype, "_Lasso.pdf", sep =""),   # The directory you want to save the file in
    width = 4.5, # The width of the plot in inches
    height = 6) # The height of the plot in inches

{
  auc_365 = ROC$AUC[[1]]
  auc_1095 = ROC$AUC[[2]]
  auc_1825 = ROC$AUC[[3]]
  dat <- data.frame(tpr365 = ROC$TP[,1],
                    fpr365 = ROC$FP[,1],
                    tpr1095 = ROC$TP[,2],
                    fpr1095 = ROC$FP[,2],
                    tpr1825 = ROC$TP[,3],
                    fpr1825 = ROC$FP[,3])
  
  ggplot() +
    geom_line(data = dat,aes(x = fpr365, y = tpr365),color = "#0066E7") +
    geom_line(data = dat,aes(x = fpr1095, y = tpr1095),color = "#FFA500")+
    geom_line(data = dat,aes(x = fpr1825, y = tpr1825),color = "#DD1717")+
    geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
    theme_bw(base_size = 20)+
    annotate("text",x = .3, y = .2, hjust = 0,
             label = paste("AUC of 1 year: ", round(auc_365,2)), color = "#0066E7", size = 15/.pt)+
    annotate("text",x = .3, y = .125, hjust = 0,
             label = paste("AUC of 3 years: ", round(auc_1095,2)), color = "#FFA500", size = 15/.pt)+
    annotate("text",x = .3, y = .05, hjust = 0,
             label = paste("AUC of 5 years: ", round(auc_1825,2)), color = "#DD1717", size = 15/.pt)+
    scale_x_continuous(name = "FPR")+
    scale_y_continuous(name = "TPR")+
    theme(axis.title.y = element_text(size = 20))
}

dev.off()

coxdata <- data # test data
colnames(coxdata)[1:2] <- c("OS.time", "OS")
#================================================================
