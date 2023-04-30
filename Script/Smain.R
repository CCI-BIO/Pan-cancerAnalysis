# Katana
# Clear the environment
rm(list = ls())

# Load necessary libraries if any
# Katana
library(TCGAbiolinks)
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
library("optparse")
library(grid)
library(futile.logger)
library(VennDiagram)
library(ezcox)
library(glmnet)
library(survival)
library(edgeR) # Loading required package: limma
# Gadi

option_list = list(
  make_option(c("-r", "--rootDir"), type = "character"),
  make_option(c("-o", "--outDir"), type = "character"),
  make_option(c("-t", "--controlDir"), type = "character"),
  make_option(c("-c", "--cancertype"), type = "character"),
  make_option(c("-p", "--primarySite"), type = "character"),
  make_option(c("-b", "--bodySite"), type = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#---------------------------------------
# Set environment variables if any
# Please remember to create necessary folders
# Katana
rootDir <- opt$rootDir
outDir <- opt$outDir
controlDir <- opt$controlDir
cancertype <- opt$cancertype
primarySite <- c(opt$primarySite)
# GBM
# LGG
# BLCA
# dir.create(file.path("/srv/scratch/z3538133/002NetworkAnalysis/Data", "PAAD"), showWarnings = FALSE)
# rootDir="/srv/scratch/z3538133/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# outDir="/srv/scratch/z3538133/002NetworkAnalysis/Data/LGG"
# controlDir="/srv/scratch/z3538133/002NetworkAnalysis/control"
# cancertype="LGG"
# primarySite=""
if (primarySite == "Brain") {
  bodySite <- "Brain - Cortex"  
} else {
  bodySite <- NA
}
if (primarySite == "Adrenal") {
  primarySite <- "Adrenal Gland"
}

#---------------------------------------

# Set working directory
setwd(rootDir)

# Include the script of functions
source(paste(rootDir, "/Script/ProposedMethod_Functions.R", sep=""))

# Create folders
# dir.create(file.path(paste(rootDir, "/Data", sep = ""), cancertype), showWarnings = FALSE)

#================================================================
# (1) Update limma, save to file DEG_limmavoom.Rdata
#================================================================
f <- paste(outDir, "/", cancertype, "_output.txt", sep = "")
write("Starting",file=f,append=FALSE)

# Load data
load(paste(outDir, '/gtex.', cancertype, '.meta.RData', sep = ""))
head(gtex.PAN.meta)
nrow(gtex.PAN.meta)
# > head(gtex.PAN.meta)
# Sample_ID Condition   Batch
# 1 GTEX.111CU.0126.SM.5GZWZ    Normal Batch 1
# 2 GTEX.111YS.0126.SM.5987T    Normal Batch 1
# 3 GTEX.1122O.0326.SM.5H124    Normal Batch 1
# 4 GTEX.11DXX.0126.SM.5EGH7    Normal Batch 1
# 5 GTEX.11DXY.1626.SM.5H12L    Normal Batch 1
# 6 GTEX.11DXZ.0226.SM.5EGGZ    Normal Batch 1
# > nrow(gtex.PAN.meta)
# [1] 202

load(paste(outDir, '/gtex.', cancertype, '.raw.counts.RData', sep = ""))
gtex.PAN.raw.counts[1:5,1:3]
nrow(gtex.PAN.raw.counts)
ncol(gtex.PAN.raw.counts)
# > gtex.PAN.raw.counts[1:5,1:3]
# gene TCGA.OR.A5K3.01 TCGA.OR.A5J2.01
# 1 5_8S_rRNA         0.00000         0.00000
# 2   5S_rRNA         0.00000         0.00000
# 3       7SK         0.00000         0.00000
# 4      A1BG        27.00087        19.99975
# 5  A1BG-AS1         2.58010        34.93954
# > nrow(gtex.PAN.raw.counts)
# [1] 58581
# > ncol(gtex.PAN.raw.counts)
# [1] 203

# 1.1) Analysis with limma

# assign gene column as rowname
rownames(gtex.PAN.raw.counts) <- gtex.PAN.raw.counts[,"gene"]
gtex.PAN.raw.counts[ ,"gene"] <- NULL

group_list <- gtex.PAN.meta[ ,2]

exprSet <- gtex.PAN.raw.counts[ ,sort(names(gtex.PAN.raw.counts))] # to move GTEx before TCGA

# dge <- DGEList(counts = exprSet, group=group_list) # Create DGEList object
dge <- DGEList(counts = exprSet) # Create DGEList object
# Calculate normalization factors
dge <- calcNormFactors(dge)
# Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
d <- dge[-drop,] 
dim(d) # number of genes left, row number
# > dim(d) # number of genes left
# [1] 23079   202
# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
group_list <- factor(group_list)
design <- model.matrix(~0 + group_list)
rownames(design) <- colnames(dge)
colnames(design) <- levels(group_list)
head(design)
# > head(design)
# Normal Tumor
# GTEX.111CU.0126.SM.5GZWZ      1     0
# GTEX.111YS.0126.SM.5987T      1     0
# GTEX.1122O.0326.SM.5H124      1     0
# GTEX.11DXX.0126.SM.5EGH7      1     0
# GTEX.11DXY.1626.SM.5H12L      1     0
# GTEX.11DXZ.0226.SM.5EGGZ      1     0

# Voom
y <- voom(d, design, plot = F)

# lmFit fits a linear model using weighted least squares for each gene
fit <- lmFit(y, design = design)

contrast.matrix <- makeContrasts("Tumor-Normal", levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

DEG_limma <- topTable(fit2, coef = 1, n = Inf)
DEG_limma <- na.omit(DEG_limma)

logFC_cutoff <- 1
DEG_limma$change <- as.factor(
  ifelse(DEG_limma$adj.P.Val < 0.05 & abs(DEG_limma$logFC) > logFC_cutoff,
         ifelse(DEG_limma$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

table(DEG_limma$change)
# > table(DEG_limma$change)
# 
# DOWN   NOT    UP
# 4476 15847  2756

f <- paste(outDir, "/", cancertype, "_output.txt", sep = "")
write("limma",file=f,append=TRUE)
write(paste(length(which(DEG_limma$change %in% c("DOWN", "UP"))), " genes' expression changed", sep = ""),file=f,append=TRUE)

DEG_limmavoom <- DEG_limma
head(DEG_limmavoom)
nrow(DEG_limmavoom)
# > head(DEG_limmavoom)
# logFC    AveExpr        t       P.Value     adj.P.Val        B
# AC009065.4  6.381303 -3.7843237 54.53280 4.165181e-125 9.612820e-121 273.0285
# RP11-40C6.2 7.256780 -2.4690511 52.35802 1.109139e-121 1.279891e-117 265.8528
# HNRNPA1P4   4.516182 -3.8366772 39.67918  6.113227e-99  4.702906e-95 214.4156
# HNRNPCP2    3.787074  0.9780526 38.62826  8.341025e-97  4.812563e-93 210.2543
# CTB-63M22.1 4.947028  2.4752223 31.55818  3.251207e-81  1.500692e-77 174.7269
# EIF5AP4     4.816664 -2.7988863 31.23088  1.949673e-80  7.499416e-77 172.6846
# change
# AC009065.4      UP
# RP11-40C6.2     UP
# HNRNPA1P4       UP
# HNRNPCP2        UP
# CTB-63M22.1     UP
# EIF5AP4         UP
# > nrow(DEG_limmavoom)
# [1] 23079

save(DEG_limmavoom, file=paste(outDir, '/DEG_limmavoom.Rdata', sep = ""))
write.csv(DEG_limmavoom,paste(outDir, '/DEG_limmavoom.csv', sep = ""), row.names = T)

#================================================================
# (2) Building the network for a specific condition
# Using STRING db
#================================================================
# Load the tumor expression data
load(paste(outDir, "/expData.RData", sep = ""))

# Get PPI network
# edges <- read_excel(paste(rootDir, "/Data/PPI.xls", sep = ""), sheet = 1)
edges <- read.csv(file = paste(rootDir, "/Data/stringInteractions.csv", sep = ""))
load(file = paste(outDir, "/DEG_limmavoom.Rdata", sep = ""))
x <- DEG_limmavoom[which(DEG_limmavoom$change %in% c("UP", "DOWN")),]

selectedGenes <- rownames(x)
# selectedGenes <- DEGlist$gene
interactions <- edges[, c(1, 2)]
colnames(interactions) <- c("cause", "effect")
interactions <- interactions[which(interactions$cause %in% colnames(expData$mRNAs)),]
interactions <- interactions[which(interactions$cause %in% selectedGenes),]
interactions <- interactions[which(interactions$effect %in% colnames(expData$mRNAs)),]
interactions <- interactions[which(interactions$effect %in% selectedGenes),]
nodes <- unique(union(interactions$cause, interactions$effect))
length(nodes)

# TFs: Download the list from http://fantom.gsc.riken.jp/5/sstar/Browse_Transcription_Factors_hg19
tfs <- read.csv(paste(rootDir, "/Data/Browse Transcription Factors hg19 - resource_browser.csv",
                      sep = ""))
i <- which(tfs$Symbol %in% nodes)
tfData <- expData$mRNAs[, tfs$Symbol[i]]

# Update cancer data of mRNAs
expData$mRNAs <- expData$mRNAs[, nodes[which(!(nodes %in% tfs$Symbol[i]))]]
mRNAsData_Cancer <-  expData$mRNAs

# Get the cancer data of miRNAs
# miRNAsData_Cancer <-  expData$miRs
miRNAsData_Cancer <-  NA

# Remove genes with over 50% samples having expression = 0
# mRNAs
n <- ncol(mRNAsData_Cancer) # number of genes
nSamples <- nrow(mRNAsData_Cancer)
l <- c()
for (i in 1:n) {
  t <- length(which(mRNAsData_Cancer[,i] == 0))
  if(t > nSamples/2) {
    l <- c(l,i)
  }
}
if(length(l) > 0) {
  mRNAsData_Cancer <- mRNAsData_Cancer[,-l]
}

# TFs
n <- ncol(tfData) # number of genes
nSamples <- nrow(tfData)
l <- c()
for (i in 1:n) {
  t <- length(which(tfData[,i] == 0))
  if(t > nSamples/2) {
    l <- c(l,i)
  }
}
if(length(l) > 0) {
  tfData <- tfData[,-l]
}

# Combine data
nomiR <- ncol(miRNAsData_Cancer) # if expData$miRs is NA, nomiR is null
nomR <- ncol(mRNAsData_Cancer)
noTF <- ncol(tfData)
if(!is.null(nomiR)) {
  cancer_data <- cbind(miRNAsData_Cancer, mRNAsData_Cancer, tfData)  
} else {
  cancer_data <- cbind(mRNAsData_Cancer, tfData)
}

# Free the memory
gc()

# Build the network
cancer_network <- buildNetwork(interactions, nomiR, nomR, noTF, cancer_data, rootDir, outDir,
                               usingpVal = TRUE, cutoff = 0.05)

# Save the network
write.csv(cancer_network, paste(outDir, "/cancer_network.csv", sep = ""), row.names = FALSE)

# # remove duplicated
# subtypes <- PanCancerAtlas_subtypes()
# subtypes <- unique(subtypes$cancer.type)
# subtypes <- c(subtypes, "PAAD")
# subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
# subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
# n <- length(subtypes)
# 
# for (i in 1:n) {
#   
#   cancertype <- subtypes[i]
#   outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
#   
#   # Network
#   network <- read.csv(paste(outDir, "/cancer_network.csv", sep = ""))
#   
#   network <- unique(network)
#   
#   # Save the network
#   write.csv(network, paste(outDir, "/cancer_network.csv", sep = ""), row.names = FALSE)
# }

# # Analyse network
# noEdge <- analyseNetwork(0,nomR, noTF, cancer_network, cancer_data,
#                paste(outDir, "/cancer_network_analysis.txt", sep = ""))
# 
# f <- paste(outDir, "/", cancertype, "_output.txt", sep = "")
# write("Number of nodes",file=f,append=TRUE)
# write(paste(nomR + noTF, " nodes", sep = ""),file=f,append=TRUE)
# write("Number of edges",file=f,append=TRUE)
# write(paste(noEdge, " edges", sep = ""),file=f,append=TRUE)
