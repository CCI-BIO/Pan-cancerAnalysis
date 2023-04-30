# Local
# Clear the environment
rm(list = ls())

# Load necessary libraries if any
require(data.table)
library(RColorBrewer)
library(edgeR)
library(survival)
library(wesanderson)
library(TCGAbiolinks)
library(ggplot2)
library(tidyverse)
library(readr)
library(reshape2)
library(maftools)
library(CancerSubtypes)
library("SNFtool")
library(ezcox)
library(g3viz)

library(datasets)
library(taRifx)
library(xtable)
library(UpSetR)
library(STRINGdb)
library(writexl)
library(dplyr)
library(mygene)
library(stringr)
library(reshape)
library(plyr)
library(scales)
library(ggpubr)
library(gridExtra)
library(openxlsx)
library(hrbrthemes)
library(viridis)
library(plotly)
library(heatmaply)
library(forcats)
library(glmnet)
library(RTCGAToolbox)
library(reticulate)
library("biomaRt")
library(EnsDb.Hsapiens.v79)
library(GGally)
library(network)
library(sna)
library(umap)

#---------------------------------------
# Set environment variables if any
# Please remember to create necessary folders
# rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"
rootDir="/srv/scratch/z3538133/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# outDir="C:/Users/vpham/Documents/002NetworkAnalysis/Data/LGG"
# cancertype="LGG"
#---------------------------------------

#---------------------------------------
genome<-"hg38" # insert genome used either mm10, hg19, or hg38
if(genome == "mm10"){
  library(org.Mm.eg.db)
  orgSpia<-"mmu"
  egSymbol<-toTable(org.Mm.egSYMBOL)
  genomeFull<-"Mus musculus"
}else if(genome == "hg19" || genome == "hg38"){
  library(org.Hs.eg.db)
  orgSpia<-"hsa"
  egSymbol<-toTable(org.Hs.egSYMBOL)
  genomeFull<-"Homo sapiens"
}
#---------------------------------------

# Set working directory
setwd(rootDir)

# Include the script of functions
source(paste(rootDir, "/Script/ProposedMethod_Functions.R", sep=""))

# qsub -I -l select=1:ncpus=2:mem=24gb,walltime=9:00:00

# module load gcc/7.5.0
# module load R/4.0.2-gcc7
# module unload gcc/7.5.0
# module load gcc/8.4.0

# qsub -I -l select=1:ncpus=2:mem=24gb,walltime=9:00:00
# module load gcc/12.2.0
# module load r/4.2.2
# module load gsl/2.7.1
# R

# #================================================================
# # (-2) Table for data
# #================================================================
# 
# subtypes <- PanCancerAtlas_subtypes()
# subtypes <- unique(subtypes$cancer.type)
# subtypes
# # > subtypes
# # [1] "ACC"  "AML"  "BLCA" "BRCA" "LGG"  "GBM"  "ESCA" "COAD" "STAD" "READ" "HNSC" "KICH" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "OVCA" "PCPG" "PRAD" "SKCM"
# # [22] "THCA" "UCEC" "UCS" 
# subtypes <- c(subtypes, "PAAD")
# subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
# subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
# subtypes
# # > subtypes
# # [1] "ACC"  "AML"  "BLCA" "BRCA" "COAD" "ESCA" "GBM"  "KICH" "KIRC" "KIRP" "LGG"  "LIHC" "LUAD" "LUSC" "OVCA" "PAAD" "PCPG" "PRAD" "SKCM" "STAD" "THCA"
# # [22] "UCEC" "UCS" 
# 
# names <- c("Adrenocortical carcinoma", "Acute myeloid leukemia", "Bladder Urothelial Carcinoma", "Breast invasive carcinoma", "Colon adenocarcinoma", "Esophageal carcinoma",
#            "Glioblastoma multiforme", "Kidney Chromophobe", "Kidney renal clear cell carcinoma", "Kidney renal papillary cell carcinoma", "Brain Lower Grade Glioma", "Liver hepatocellular carcinoma",
#            "Lung adenocarcinoma", "Lung squamous cell carcinoma", "Serous ovarian carcinoma", "Pancreatic Cancer", "Pheochromocytoma and Paraganglioma", "Prostate adenocarcinoma",
#            "Skin Cutaneous Melanoma", "Stomach adenocarcinoma", "Thyroid carcinoma", "Uterine Corpus Endometrial Carcinoma", "Uterine Carcinosarcoma")
# tumour <- c(76, 172, 144, 1212, 279, 179,
#             166, 91, 508, 189, 348, 239,
#             274, 192, 418, 154, 179, 373,
#             331, 409, 560, 182, 57)
# normal <- c(126, 444, 9, 178, 307, 652,
#             105, 28, 28, 28, 105, 110,
#             288, 288, 88, 167, 126, 100,
#             811, 174, 279, 78, 78)
# 
# n <- length(subtypes)
# t <- matrix(NA, nrow = n, ncol = 4)
# colnames(t) <- c('Code of cancer type', 'Cancer type', 'TCGA primary tumor', 'GTEx normal tissue')
# for (i in 1:n) {
#   cancertype <- subtypes[i]
#   t[i,1] <- cancertype
#   t[i,2] <- names[i]
#   t[i,3] <- as.numeric(tumour[i])
#   t[i,4] <- as.numeric(normal[i])
# }
# 
# t[,3]<-format(as.numeric(t[,3]),big.mark=",")
# t[,4]<-format(as.numeric(t[,4]),big.mark=",")
# 
# outDir <- paste(rootDir, "/Data/Output/Table", sep = "")
# write.xlsx(as.data.frame(t), paste(outDir, "/DataTable.xlsx", sep = ""))
# 
# #================================================================

#================================================================
# (-1) Latex table for differentially expressed genes
#================================================================

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

n <- length(subtypes)
t <- matrix(NA, nrow = n, ncol = 2)
colnames(t) <- c('Cancer type', 'Number of differentially expressed genes')
for (i in 1:n) {
  cancertype <- subtypes[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(file=paste(outDir, '/DEG_limmavoom.Rdata', sep = "")) # return DEG_limmavoom
  numGenes <- length(which(DEG_limmavoom$change %in% c("DOWN", "UP")))
  t[i,1] <- cancertype
  t[i,2] <- as.numeric(numGenes)
}

t[,2]<-format(as.numeric(t[,2]),big.mark=",")
print(latex.table.by(t, format.args = list(digits = 2, format = c("s","d"))), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

#================================================================

#================================================================
# (0) Latex table for networks
#================================================================

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes
# > subtypes
# [1] "ACC"  "AML"  "BLCA" "BRCA" "LGG"  "GBM"  "ESCA" "COAD" "STAD" "READ" "HNSC" "KICH" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "OVCA" "PCPG" "PRAD" "SKCM"
# [22] "THCA" "UCEC" "UCS" 
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
subtypes
# > subtypes
# [1] "ACC"  "AML"  "BLCA" "BRCA" "COAD" "ESCA" "GBM"  "KICH" "KIRC" "KIRP" "LGG"  "LIHC" "LUAD" "LUSC" "OVCA" "PAAD" "PCPG" "PRAD" "SKCM" "STAD" "THCA"
# [22] "UCEC" "UCS" 

n <- length(subtypes)
t <- matrix(NA, nrow = n, ncol = 4)
colnames(t) <- c('Cancer type', 'Number of nodes', 'Number of edges', 'Number of critical nodes')
for (i in 1:n) {
  cancertype <- subtypes[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  # net <- read.table(paste(outDir, "/cancer_network_analysis.txt", sep = ""), header = FALSE, sep = "\t", dec = ".", fill = TRUE)
  net <- read.csv(paste(outDir, "/cancer_network.csv", sep = ""))
  t[i,1] <- cancertype
  t[i,2] <- length(unique(c(net[,1], net[,2])))
  t[i,3] <- nrow(net)
  cri <- read.csv(paste(outDir, "/critical_nodes.csv", sep = ""), header = TRUE, sep = ",", dec = ".")
  t[i,4] <- nrow(cri)
}

t[,2]<-format(as.numeric(t[,2]),big.mark=",")
t[,3]<-format(as.numeric(t[,3]),big.mark=",")
t[,4]<-format(as.numeric(t[,4]),big.mark=",")
print(latex.table.by(t, format.args = list(digits = 2, format = c("s","d","d","d"))), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

#================================================================

#================================================================
# (1) Latex table for critical copper genes
#================================================================

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes
# > subtypes
# [1] "ACC"  "AML"  "BLCA" "BRCA" "LGG"  "GBM"  "ESCA" "COAD" "STAD" "READ" "HNSC" "KICH" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "OVCA" "PCPG" "PRAD" "SKCM"
# [22] "THCA" "UCEC" "UCS" 
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
subtypes
# > subtypes
# [1] "ACC"  "AML"  "BLCA" "BRCA" "COAD" "ESCA" "GBM"  "KICH" "KIRC" "KIRP" "LGG"  "LIHC" "LUAD" "LUSC" "OVCA" "PAAD" "PCPG" "PRAD" "SKCM" "STAD" "THCA"
# [22] "UCEC" "UCS" 

n <- length(subtypes)
t <- matrix(NA, nrow = n, ncol = 3)
colnames(t) <- c('Cancer type', 'Number of CCGs', 'Critical cuproplasia-related genes (CCGs)')
for (i in 1:n) {
  cancertype <- subtypes[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(paste(outDir, "/unicox.Rdata", sep = ""))
  t[i,1] <- cancertype
  t[i,2] <- nrow(unicox)
  genes <- unicox$Variable
  # genes <- genes[order(genes, decreasing = FALSE)]
  t[i,3] <- paste(genes, collapse=", ")
}

print(latex.table.by(t), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)
#   then add \usepackage{multirow} to the preamble of your LaTeX document
#   for longtable support, add ,tabular.environment='longtable' to the print command (plus add in ,floating=FALSE), then \usepackage{longtable} to the LaTeX preamble

#================================================================

#================================================================
# (2) Check common genes in the results
#================================================================

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes
# > subtypes
# [1] "ACC"  "AML"  "BLCA" "BRCA" "LGG"  "GBM"  "ESCA" "COAD" "STAD" "READ" "HNSC" "KICH" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "OVCA" "PCPG" "PRAD" "SKCM"
# [22] "THCA" "UCEC" "UCS" 
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

n <- length(subtypes)
for (i in 1:n) {
  cancertype <- subtypes[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(paste(outDir, "/unicox.Rdata", sep = ""))
  genes <- unicox$Variable
  if(i == 1) {
    allgenes <- genes
  } else {
    allgenes <- c(allgenes, genes)  
  }
}

allgenes <- unique(allgenes)
nr <- length(allgenes)
nc <- length(subtypes)
dat <- matrix(0, nrow = nr, ncol = nc)
row.names(dat) <- allgenes
colnames(dat) <- subtypes
for (i in 1:nc) {
  cancertype <- subtypes[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(paste(outDir, "/unicox.Rdata", sep = ""))
  for (j in 1:nr) {
    if(row.names(dat)[j] %in% unicox$Variable) {
      dat[j,i] <- 1
    }
  }
}

criticalCopperGenes <- as.data.frame(dat)
rownames(criticalCopperGenes)
nr # number of critical copper genes
nc # number of cancer types
# > rownames(criticalCopperGenes)
# [1] "CDK1"     "COX17"    "DBH"      "SLC11A2"  "MAP1LC3A" "ALB"
# [7] "ANG"      "ANKRD9"   "AP1S1"    "ARF1"     "APC"      "TMPRSS6"
# [13] "GPC1"     "LOXL1"    "ATOX1"    "MAPT"     "CASP3"    "XIAP"
# [19] "GSK3B"    "AOC3"     "S100A12"  "FOXO1"    "JUN"      "SNCA"
# [25] "SORD"     "APP"      "PRNP"     "ATP6AP1"  "F5"       "ADAM10"
# [31] "ADAM17"   "COA6"     "MT-CO2"   "STEAP3"   "CYP1A1"   "CP"
# [37] "MMGT1"    "BACE1"    "AP1B1"    "TP53"     "ATP7A"    "SP1"
# [43] "CCND1"    "MT-CO1"   "PRND"     "AANAT"    "AQP1"     "MT1X"
# [49] "IL1A"     "SPATA5"   "COMMD1"   "F8"       "SUMF1"    "XAF1"
# [55] "HEPH"     "PARK7"    "LCAT"
# > nr # number of critical copper genes
# [1] 57
# > nc # number of cancer types
# [1] 23

x <- criticalCopperGenes
x$frequency <- rowSums(criticalCopperGenes)
x <- x[order(x$frequency, decreasing = TRUE),]

print(paste(row.names(x), collapse=", "))
# [1] "CDK1, ALB, AP1S1, CASP3, MAP1LC3A, SNCA, TMPRSS6, MAPT, GSK3B, JUN, APP, CYP1A1, COX17, XIAP, TP53, FOXO1, COA6, ARF1, GPC1, AOC3, SORD, PRNP, F5, ATP7A, SP1, MT-CO1, DBH, SLC11A2, ANG, S100A12, ATP6AP1, ADAM10, MT-CO2, CP, BACE1, PRND, AQP1, MT1X, IL1A, XAF1, ANKRD9, APC, LOXL1, ATOX1, ADAM17, STEAP3, MMGT1, AP1B1, CCND1, AANAT, SPATA5, COMMD1, F8, SUMF1, HEPH, PARK7, LCAT"

important.copper.genelist <- c('SLC31A1', 'MT2A', 'ATP7A', 'ATP7B', 'MT1X', 'ATOX1', 'COX17', 'MTH1', 'SOD1', 'GSS')
intersect(important.copper.genelist, row.names(x))
# > intersect(important.copper.genelist, row.names(x))
# [1] "ATP7A" "MT1X"  "ATOX1" "COX17"

# Only get critical copper genes in more than 1 cancer type
selectedCriticalCopperGenes <- x[which(x$frequency > 1),]
print(paste(row.names(selectedCriticalCopperGenes), collapse=", "))
nrow(selectedCriticalCopperGenes)
# > print(paste(row.names(selectedCriticalCopperGenes), collapse=", "))
# [1] "CDK1, ALB, AP1S1, CASP3, MAP1LC3A, SNCA, TMPRSS6, MAPT, GSK3B, JUN, APP, CYP1A1, COX17, XIAP, TP53, FOXO1, COA6, ARF1, GPC1, AOC3, SORD, PRNP, F5, ATP7A, SP1, MT-CO1, DBH, SLC11A2, ANG, S100A12, ATP6AP1, ADAM10, MT-CO2, CP, BACE1, PRND, AQP1, MT1X, IL1A, XAF1"
# > nrow(selectedCriticalCopperGenes)
# [1] 40

# Get gene activities
x$gene <- row.names(x)
m<-match(x$gene,egSymbol$symbol)
x$entrez_id<-egSymbol$gene_id[m]
#Add mygene summaries
geneSummaries<-as.data.frame(getGenes(geneids=c(x$entrez_id),fields="all",return.as="DataFrame")[,c("query","summary")])
m2<-match(x$entrez_id,geneSummaries$query)
x$summary<-geneSummaries$summary[m2]
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- x[,-which(colnames(x) == "entrez_id")]
x <- x %>% dplyr::select(frequency, everything())
x <- x %>% dplyr::select(gene, everything())
write.csv(x,fileName,row.names=F)

# Only get interested genes
selectedCriticalCopperGenes <- x[which(x$frequency > 1),]
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/topCriticalCoperGenes.csv",sep="")
write.csv(selectedCriticalCopperGenes,fileName,row.names=F)
selectedCriticalCopperGenes <- selectedCriticalCopperGenes[, c("gene", "frequency", "summary")]
selectedCriticalCopperGenes$frequency <- as.character(selectedCriticalCopperGenes$frequency)
print(latex.table.by(selectedCriticalCopperGenes), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

#================================================================

#================================================================
# (3) Get STRING db data
# local
#================================================================

# Get the data
# 9606 for Human
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=700, network_type="full", input_directory="")
string_proteins <- string_db$get_proteins()
string_proteins[1:5,1:3]
nrow(string_proteins)
# > string_proteins[1:5,1:3]
# protein_external_id preferred_name protein_size
# 1 9606.ENSP00000000233           ARF5          180
# 2 9606.ENSP00000000412           M6PR          277
# 3 9606.ENSP00000001008          FKBP4          459
# 4 9606.ENSP00000001146        CYP26B1          512
# 5 9606.ENSP00000002125        NDUFAF7          441
# > nrow(string_proteins)
# [1] 19566
stringInteractions <- string_db$get_interactions(string_proteins$protein_external_id)

# Convert to gene names
t <- merge(stringInteractions, string_proteins, by.x = "from", by.y = "protein_external_id")
t <- t[,c(2,3,4)]
colnames(t)[3] <- "from"

stringInteractions <- merge(t, string_proteins, by.x = "to", by.y = "protein_external_id")
stringInteractions <- stringInteractions[,c(2,3,4)]
colnames(stringInteractions)[3] <- "to"
stringInteractions <- stringInteractions %>% dplyr::select(to, everything())
stringInteractions <- stringInteractions %>% dplyr::select(from, everything())
stringInteractions <- stringInteractions[,c(1,2)]

head(stringInteractions)
nrow(stringInteractions)
# > head(stringInteractions)
# from      to
# 1 FKBP4   PPP5C
# 2 FKBP4   PPP5C
# 3  RALA  RALBP1
# 4  RALA  RALBP1
# 5 XYLT2 B4GALT7
# 6 XYLT2 B4GALT7
# > nrow(stringInteractions)
# [1] 505968

# Duplicated, get unique
stringInteractions <- unique(stringInteractions)
head(stringInteractions)
nrow(stringInteractions)
# > head(stringInteractions)
# from        to
# 1   FKBP4     PPP5C
# 3    RALA    RALBP1
# 5   XYLT2   B4GALT7
# 7  RB1CC1 GABARAPL2
# 9    CFTR     PSMA4
# 11    CD4      LCP2
# > nrow(stringInteractions)
# [1] 252956

# Save file
write.csv(stringInteractions,paste(rootDir, "/Data/stringInteractions.csv", sep = ""), row.names = FALSE)
write.csv(string_proteins[,1:3],paste(rootDir, "/Data/string_proteins.csv", sep = ""), row.names = FALSE)

#================================================================

#================================================================
# (4) BRCA - RNA Seq data
#================================================================

outDir="C:/Users/vpham/Documents/002NetworkAnalysis/Data/BRCA"

# TEPA8
fileName <- paste(outDir, "/RNAseq/DrugvsControl.MB231.TEPA8.txt", sep = "")
deg_expTab8 <- read.table(fileName, sep = "\t", header = TRUE)
deg_expTab8 <- deg_expTab8[,c(1,2,7,5,6)]
deg_expTab8 <- deg_expTab8[order(abs(deg_expTab8$PValue), decreasing = FALSE),]
print(latex.table.by(deg_expTab8, digits = c(0,0,3,3,-3,-3)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

# TEPA24
fileName <- paste(outDir, "/RNAseq/DrugvsControl.MB231.TEPA24.txt", sep = "")
deg_expTab24 <- read.table(fileName, sep = "\t", header = TRUE)
deg_expTab24 <- deg_expTab24[,c(1,2,7,5,6)]
deg_expTab24 <- deg_expTab24[order(abs(deg_expTab24$PValue), decreasing = FALSE),]
print(latex.table.by(deg_expTab24, digits = c(0,0,3,3,-3,-3)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

#================================================================

#================================================================
# (5) Compare gene expression
#================================================================

# Run on server

# Link the files
# Katana
# cancertype="ACC"
# ln -s /srv/scratch/z3538133/001pancancer/pan/data/$cancertype/gtex.$cancertype.raw.counts.RData /srv/scratch/z3538133/002NetworkAnalysis/Data/$cancertype/gtex.$cancertype.raw.counts.RData

# Create a plot
prepareBoxPlotData=function(gene) {
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- unique(subtypes$cancer.type)
  subtypes <- c(subtypes, "PAAD")
  subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
  subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
  n <- length(subtypes)
  Cancer_type <- c()
  Type <- c()
  Expression <- c()
  for (i in 1:n) {
    cancertype <- subtypes[i]
    outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
    load(paste(outDir, "/gtex.", cancertype, ".raw.counts.RData", sep = "")) # return gtex.PAN.raw.counts
    # transform
    dat <- gtex.PAN.raw.counts
    rownames(dat) <- dat$gene
    dat <- dat[,-1]
    dat <- dat + 1
    dat <- log(dat, base = 2)
    iGene <- which(rownames(dat) == gene)
    r <- str_detect(colnames(dat), "TCGA")
    nTumor <- length(which(r == TRUE))
    nNormal <- ncol(dat) - nTumor
    
    # Add * if being critical copper gene
    outDir <- paste(rootDir, "/Data/Output", sep = "")
    fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
    selectedCriticalCopperGenes <- read.csv(fileName)
    if(selectedCriticalCopperGenes[which(selectedCriticalCopperGenes$gene == gene),cancertype] == 1){
      cancertype <- paste("*", cancertype, sep = "")
    }
    
    Cancer_type <- c(Cancer_type, rep(cancertype, nTumor + nNormal))
    Type <- c(Type, rep(c("Tumour", "Normal"), times = c(nTumor, nNormal)))
    Expression <- c(Expression, as.numeric(dat[iGene,]))
  }
  data=data.frame(Cancer_type, Type ,  Expression)
  
  return(data)
}

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/topCriticalCoperGenes.csv",sep="")
selectedCriticalCopperGenes <- read.csv(fileName)
selectedCriticalCopperGenes <- selectedCriticalCopperGenes[1:10,]
n <- nrow(selectedCriticalCopperGenes)
dir.create(file.path(outDir, "boxplot"), showWarnings = FALSE)

selectedCriticalCopperGenes$gene
# > selectedCriticalCopperGenes$gene
# [1] "CDK1"     "ALB"      "AP1S1"    "CASP3"    "MAP1LC3A" "SNCA"
# [7] "TMPRSS6"  "MAPT"     "GSK3B"    "JUN"
n
# > n
# [1] 10

# i <- 1

# Run by hand
f <- paste(outDir, "/boxplot/", "combine.png", sep = "")
png(file = f, width = 595, height =  842)

i <- 1
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g1 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 2
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g2 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 3
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g3 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 4
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g4 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 5
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g5 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 6
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g6 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 7
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g7 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 8
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g8 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 9
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g9 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

i <- 10
gene <- selectedCriticalCopperGenes$gene[i]
data <- prepareBoxPlotData(gene)
g10 <- ggplot(data, aes(x=Cancer_type, y=Expression, color=Type)) + 
  geom_boxplot(lwd=0.25, fatten = 0.25, outlier.size = 0.25) +
  scale_color_manual(values=c("#0066E7", "#DD1717")) + 
  labs(y= paste(gene, " expression", sep = ""), x = "") + 
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, 
          ncol = 2, nrow = 5, legend = "bottom", common.legend = TRUE)

dev.off()

#================================================================

#================================================================
# (6) Pathway
# Local
#================================================================

rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# Set working directory
setwd(rootDir)
# Include the script of functions
source(paste(rootDir, "/Script/ProposedMethod_Functions.R", sep=""))

# Read critical copper genes
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# All 57 critical copper genes
geneList57 <- x$gene

# 40 popular critical cooper genes
geneList40 <- x[which(x$frequency > 1),]$gene

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# > head(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# > nrow(uni_gene)
# [1] 35

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]
print("35 genes")
print(paste(uni_gene$x, collapse = ", "))
print("18 genes")
print(paste(geneList18, collapse = ", "))
print("12 genes")
print(paste(geneList12, collapse = ", "))
# > print("35 genes")
# [1] "35 genes"
# > print(paste(uni_gene$x, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, MAP1LC3A, SNCA, TMPRSS6, MAPT, GSK3B, JUN, APP, CYP1A1, COX17, XIAP, FOXO1, ARF1, GPC1, AOC3, SORD, PRNP, ATP7A, SP1, MT-CO1, DBH, SLC11A2, ANG, S100A12, ATP6AP1, ADAM10, MT-CO2, CP, PRND, AQP1, MT1X, IL1A, XAF1"
# > print("18 genes")
# [1] "18 genes"
# > print(paste(geneList18, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, ARF1, GPC1, SORD, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"
# > print("12 genes")
# [1] "12 genes"
# > print(paste(geneList12, collapse = ", "))
# [1] "MAP1LC3A, SNCA, MAPT, JUN, CYP1A1, AOC3, PRNP, DBH, S100A12, AQP1, MT1X, XAF1"

geneList30 <- c(geneList18, geneList12)
cat(geneList30, sep = '\n')
# > cat(geneList30, sep = '\n')
# CDK1
# AP1S1
# CASP3
# TMPRSS6
# GSK3B
# APP
# COX17
# XIAP
# ARF1
# GPC1
# SORD
# ATP7A
# SP1
# MT-CO1
# SLC11A2
# ATP6AP1
# ADAM10
# CP
# MAP1LC3A
# SNCA
# MAPT
# JUN
# CYP1A1
# AOC3
# PRNP
# DBH
# S100A12
# AQP1
# MT1X
# XAF1

# Do enrichment analysis

# Analyse the results

# Biological process
GO_process <- read.table(paste(outDir, "/Enrich/GO_Biological_Process_2021_table.txt", sep=""),
                         as.is = TRUE, sep = "\t", header = TRUE, quote="")
GO_process <- GO_process[, c(1,2,4,9)]
colnames(GO_process) <- c("Term", "Overlap", "Adjusted p-value", "Genes")
GO_process <- GO_process[which(GO_process[,3] < 0.05),]
GO_process <- GO_process[order(GO_process[,3], decreasing = FALSE),]
write.csv(GO_process, file = paste(outDir, "/Enrich/GO_Biological_Process_2021_table_Cutoff.csv",
                                   sep=""), row.names = FALSE, quote=TRUE)
r <- prepareDataForClustergram(GO_process, termTop = 20, geneTop  = 0)
rr <- t(r)
rr <- cbind(rownames(rr),rr)
colnames(rr) <- rr[1,]
rr <- rr[-1,]
colnames(rr)[1] <- "Term"
write.csv(rr, file = paste(outDir, "/Enrich/GO_Biological_Process_2021_table_Heatmap.csv",
                          sep=""), row.names = FALSE, quote=TRUE)

# Output
# Image
dat <- read.csv(file = paste(outDir, "/Enrich/GO_Biological_Process_2021_table_Heatmap.csv", sep=""))
f <- paste(outDir, "/Enrich/GO_Biological_Process_2021.pdf", sep = "")
pdf(file = f,  width = 12, height = 8, onefile=FALSE)
l <- dat$Term
drawClustergram2(dat, "", l)
dev.off()
# Table
GO_process[,4] <- gsub("[;]","; ",GO_process[,4])
print(latex.table.by(GO_process[1:20,], digits = c(0,0,0,-3,0)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

# Molecular function
GO_process <- read.table(paste(outDir, "/Enrich/GO_Molecular_Function_2021_table.txt", sep=""),
                         as.is = TRUE, sep = "\t", header = TRUE, quote="")
GO_process <- GO_process[, c(1,2,4,9)]
colnames(GO_process) <- c("Term", "Overlap", "Adjusted p-value", "Genes")
GO_process <- GO_process[which(GO_process[,3] < 0.05),]
GO_process <- GO_process[order(GO_process[,3], decreasing = FALSE),]
write.csv(GO_process, file = paste(outDir, "/Enrich/GO_Molecular_Function_2021_table_Cutoff.csv",
                                   sep=""), row.names = FALSE, quote=TRUE)
r <- prepareDataForClustergram(GO_process, termTop = 20, geneTop  = 0)
rr <- t(r)
rr <- cbind(rownames(rr),rr)
colnames(rr) <- rr[1,]
rr <- rr[-1,]
colnames(rr)[1] <- "Term"
write.csv(rr, file = paste(outDir, "/Enrich/GO_Molecular_Function_2021_table_Heatmap.csv",
                           sep=""), row.names = FALSE, quote=TRUE)

# Output
# Image
dat <- read.csv(file = paste(outDir, "/Enrich/GO_Molecular_Function_2021_table_Heatmap.csv", sep=""))
f <- paste(outDir, "/Enrich/GO_Molecular_Function_2021.pdf", sep = "")
pdf(file = f,  width = 12, height = 8, onefile=FALSE)
l <- dat$Term
drawClustergram2(dat, "", l)
dev.off()
# Table
GO_process[,4] <- gsub("[;]","; ",GO_process[,4])
print(latex.table.by(GO_process[1:20,], digits = c(0,0,0,-3,0)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

# KEGG
top <- 20
GO_process <- read.table(paste(outDir, "/Enrich/KEGG_2021_Human_table.txt", sep=""),
                         as.is = TRUE, sep = "\t", header = TRUE, quote="")
GO_process <- GO_process[, c(1,2,4,9)]
colnames(GO_process) <- c("Term", "Overlap", "Adjusted p-value", "Genes")
GO_process <- GO_process[which(GO_process[,3] < 0.05),]
GO_process <- GO_process[order(GO_process[,3], decreasing = FALSE),]
write.csv(GO_process, file = paste(outDir, "/Enrich/KEGG_2021_Human_table_Cutoff.csv",
                                   sep=""), row.names = FALSE, quote=TRUE)
r <- prepareDataForClustergram(GO_process, termTop = top, geneTop  = 0)
rr <- t(r)
rr <- cbind(rownames(rr),rr)
colnames(rr) <- rr[1,]
rr <- rr[-1,]
colnames(rr)[1] <- "Term"
rr[1:top,1] <- GO_process$Term[1:top]
write.csv(rr, file = paste(outDir, "/Enrich/KEGG_2021_Human_table.csv",
                           sep=""), row.names = FALSE, quote=TRUE)

# Output
# Image
top <- 20
dat <- read.csv(file = paste(outDir, "/Enrich/KEGG_2021_Human_table.csv", sep=""))
f <- paste(outDir, "/Enrich/KEGG_2021_Human_table.pdf", sep = "")
pdf(file = f,  width = 12, height = 8, onefile=FALSE)
# > dat$Term
# [1] "Ferroptosis"                                                "Epithelial cell signaling in Helicobacter pylori infection"
# [3] "Alzheimer disease"                                          "Pathways of neurodegeneration"                             
# [5] "Mineral absorption"                                         "Human immunodeficiency virus 1 infection"                  
# [7] "Lipid and atherosclerosis"                                  "Colorectal cancer"                                         
# [9] "IL-17 signaling pathway"                                    "Lysosome"                                                  
# [11] "Pathways in cancer"                                         "Measles"                                                   
# [13] "Apoptosis"                                                  "Tyrosine metabolism"                                       
# [15] "Breast cancer"                                              "Non-alcoholic fatty liver disease"                         
# [17] "Hepatitis B"                                                "Vibrio cholerae infection"                                 
# [19] "Kaposi sarcoma-associated herpesvirus infection"            "Pathogenic Escherichia coli infection"    

dat$Term[2] <- "Epithelial cell signaling"
dat$Term[19] <- "Herpesvirus infection"
dat$Term
# > dat$Term
# [1] "Ferroptosis"                              "Epithelial cell signaling"                "Alzheimer disease"                       
# [4] "Pathways of neurodegeneration"            "Mineral absorption"                       "Human immunodeficiency virus 1 infection"
# [7] "Lipid and atherosclerosis"                "Colorectal cancer"                        "IL-17 signaling pathway"                 
# [10] "Lysosome"                                 "Pathways in cancer"                       "Measles"                                 
# [13] "Apoptosis"                                "Tyrosine metabolism"                      "Breast cancer"                           
# [16] "Non-alcoholic fatty liver disease"        "Hepatitis B"                              "Vibrio cholerae infection"               
# [19] "Herpesvirus infection"                    "Pathogenic Escherichia coli infection"   

l <- dat$Term
drawClustergram2(dat, "", l)
dev.off()
# Table
GO_process[,4] <- gsub("[;]","; ",GO_process[,4])
print(latex.table.by(GO_process[1:20,], digits = c(0,0,0,-3,0)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

# Reactome
top <- 20
GO_process <- read.table(paste(outDir, "/Enrich/Reactome_2022_table.txt", sep=""),
                         as.is = TRUE, sep = "\t", header = TRUE, quote="")
GO_process <- GO_process[, c(1,2,4,9)]
colnames(GO_process) <- c("Term", "Overlap", "Adjusted p-value", "Genes")
GO_process <- GO_process[which(GO_process[,3] < 0.05),]
GO_process <- GO_process[order(GO_process[,3], decreasing = FALSE),]
write.csv(GO_process, file = paste(outDir, "/Enrich/Reactome_2022_table_Cutoff.csv",
                                   sep=""), row.names = FALSE, quote=TRUE)
r <- prepareDataForClustergram(GO_process, termTop = top, geneTop  = 0)
rr <- t(r)
rr <- cbind(rownames(rr),rr)
colnames(rr) <- rr[1,]
rr <- rr[-1,]
colnames(rr)[1] <- "Term"
rr[1:top,1] <- GO_process$Term[1:top]
write.csv(rr, file = paste(outDir, "/Enrich/Reactome_2022_table.csv",
                           sep=""), row.names = FALSE, quote=TRUE)

# Output
# Image
dat <- read.csv(file = paste(outDir, "/Enrich/Reactome_2022_table.csv", sep=""))
f <- paste(outDir, "/Enrich/Reactome_2022_table.pdf", sep = "")
pdf(file = f,  width = 12, height = 8, onefile=FALSE)
dat$Term
# > dat$Term
# [1] "Activation Of Caspases Thru Apoptosome-Mediated Cleavage R-HSA-111459"                          
# [2] "SMAC (DIABLO) Binds To IAPs R-HSA-111463"                                                       
# [3] "SMAC, XIAP-regulated Apoptotic Response R-HSA-111469"                                           
# [4] "Iron Uptake And Transport R-HSA-917937"                                                         
# [5] "Caspase-mediated Cleavage Of Cytoskeletal Proteins R-HSA-264870"                                
# [6] "Advanced Glycosylation Endproduct Receptor Signaling R-HSA-879415"                              
# [7] "Cytochrome C-Mediated Apoptotic Response R-HSA-111461"                                          
# [8] "Amyloid Fiber Formation R-HSA-977225"                                                           
# [9] "MyD88 Cascade Initiated On Plasma Membrane R-HSA-975871"                                        
# [10] "Apoptotic Factor-Mediated Response R-HSA-111471"                                                
# [11] "Deregulated CDK5 Triggers Neurodegenerative Pathways In Alzheimers Disease Models R-HSA-8862803"
# [12] "TRAF6 Mediated Induction Of NFkB And MAP Kinases Upon TLR7/8 Or 9 Activation R-HSA-975138"      
# [13] "MyD88 Dependent Cascade Initiated On Endosome R-HSA-975155"                                     
# [14] "Toll Like Receptor 3 (TLR3) Cascade R-HSA-168164"                                               
# [15] "Toll Like Receptor 7/8 (TLR7/8) Cascade R-HSA-168181"                                           
# [16] "Insertion Of Tail-Anchored Proteins Into Endoplasmic Reticulum Membrane R-HSA-9609523"          
# [17] "Post-translational Protein Phosphorylation R-HSA-8957275"                                       
# [18] "Toll Like Receptor 9 (TLR9) Cascade R-HSA-168138"                                               
# [19] "MyD88-independent TLR4 Cascade R-HSA-166166"                                                    
# [20] "Degradation Of Extracellular Matrix R-HSA-1474228"     
for (i in 1:top) {
  strL <- nchar(dat$Term[i])
  dat$Term[i] <- str_trim(substr(dat$Term[i],strL-12,strL))
}
dat$Term
# > dat$Term
# [1] "R-HSA-111459"  "R-HSA-111463"  "R-HSA-111469"  "R-HSA-917937"  "R-HSA-264870"  "R-HSA-879415"  "R-HSA-111461"  "R-HSA-977225"  "R-HSA-975871" 
# [10] "R-HSA-111471"  "R-HSA-8862803" "R-HSA-975138"  "R-HSA-975155"  "R-HSA-168164"  "R-HSA-168181"  "R-HSA-9609523" "R-HSA-8957275" "R-HSA-168138" 
# [19] "R-HSA-166166"  "R-HSA-1474228"
l <- dat$Term
drawClustergram2(dat, "", l)
dev.off()
# Table
GO_process[,4] <- gsub("[;]","; ",GO_process[,4])
print(latex.table.by(GO_process[1:20,], digits = c(0,0,0,-3,0)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

#================================================================
# (7) Compare gene expression, identify critical copper genes up and down
#================================================================

# Run on server

# Link the files
# Katana
# cancertype="ACC"
# ln -s /srv/scratch/z3538133/001pancancer/pan/data/$cancertype/gtex.$cancertype.raw.counts.RData /srv/scratch/z3538133/002NetworkAnalysis/Data/$cancertype/gtex.$cancertype.raw.counts.RData

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
criticalCopperGenes <- read.csv(fileName)
n <- nrow(criticalCopperGenes)
criticalCopperGenes[1:5,1:6]
n
# > criticalCopperGenes[1:5,1:6]
# gene frequency ACC AML BLCA BRCA
# 1     CDK1        18   1   1    1    1
# 2      ALB        13   0   1    0    1
# 3    AP1S1        11   0   1    0    0
# 4    CASP3         9   0   0    0    1
# 5 MAP1LC3A         8   0   1    0    1
# > n
# [1] 57

r <- criticalCopperGenes[,-ncol(criticalCopperGenes)] # remove summary

for (i in 3:(ncol(criticalCopperGenes)-1)) {
  cancertype <- colnames(criticalCopperGenes)[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(paste(outDir, "/DEG_limmavoom.Rdata", sep = "")) # return DEG_limmavoom
  idx <- which(rownames(DEG_limmavoom) %in% criticalCopperGenes$gene)
  for (j in 1:length(idx)) {
    rowId <- idx[j]
    gene <- rownames(DEG_limmavoom)[rowId]
    if(r[which(r$gene == gene),i] == 1) {
      r[which(r$gene == gene),i] <- as.character(DEG_limmavoom[rowId,7])
    }
  }
}

# check r
head(r)
# > head(r)
# gene frequency ACC  AML BLCA BRCA COAD ESCA  GBM KICH KIRC KIRP  LGG LIHC
# 1     CDK1        18  UP   UP   UP   UP   UP   UP    0   UP    0    0   UP   UP
# 2      ALB        13   0 DOWN    0 DOWN DOWN    0 DOWN DOWN    0 DOWN DOWN    0
# 3    AP1S1        11   0 DOWN    0    0   UP    0    0    0    0   UP    0   UP
# 4    CASP3         9   0    0    0   UP   UP    0   UP   UP    0   UP   UP    0
# 5 MAP1LC3A         8   0 DOWN    0 DOWN    0 DOWN    0    0    0    0 DOWN    0
# 6     SNCA         8   0    0    0 DOWN DOWN    0 DOWN    0    0    0 DOWN    0
# LUAD LUSC OVCA PAAD PCPG PRAD SKCM STAD THCA UCEC UCS
# 1   UP    0   UP   UP   UP   UP   UP   UP   UP   UP   0
# 2    0    0 DOWN DOWN    0 DOWN DOWN DOWN DOWN    0   0
# 3   UP    0   UP    0   UP   UP    0   UP    0   UP  UP
# 4    0    0   UP   UP    0    0    0    0    0   UP   0
# 5    0 DOWN    0    0   UP    0 DOWN    0 DOWN    0   0
# 6    0    0 DOWN    0   UP    0   UP    0    0 DOWN   0
r$numUp <- rowSums(r == "UP")
r$numDown <- rowSums(r == "DOWN")
r[1:10,]
# > r[1:10,]
# gene frequency ACC  AML BLCA BRCA COAD ESCA  GBM KICH KIRC KIRP  LGG
# 1      CDK1        18  UP   UP   UP   UP   UP   UP    0   UP    0    0   UP
# 2       ALB        13   0 DOWN    0 DOWN DOWN    0 DOWN DOWN    0 DOWN DOWN
# 3     AP1S1        11   0 DOWN    0    0   UP    0    0    0    0   UP    0
# 4     CASP3         9   0    0    0   UP   UP    0   UP   UP    0   UP   UP
# 5  MAP1LC3A         8   0 DOWN    0 DOWN    0 DOWN    0    0    0    0 DOWN
# 6      SNCA         8   0    0    0 DOWN DOWN    0 DOWN    0    0    0 DOWN
# 7   TMPRSS6         7   0 DOWN    0    0    0    0    0 DOWN    0   UP    0
# 8      MAPT         7   0    0 DOWN    0    0 DOWN DOWN    0    0    0    0
# 9     GSK3B         7   0    0    0   UP    0   UP    0    0   UP   UP    0
# 10      JUN         6   0    0    0 DOWN    0    0   UP    0   UP    0    0
# LIHC LUAD LUSC OVCA PAAD PCPG PRAD SKCM STAD THCA UCEC UCS numUp numDown
# 1    UP   UP    0   UP   UP   UP   UP   UP   UP   UP   UP   0    18       0
# 2     0    0    0 DOWN DOWN    0 DOWN DOWN DOWN DOWN    0   0     0      13
# 3    UP   UP    0   UP    0   UP   UP    0   UP    0   UP  UP    10       1
# 4     0    0    0   UP   UP    0    0    0    0    0   UP   0     9       0
# 5     0    0 DOWN    0    0   UP    0 DOWN    0 DOWN    0   0     1       7
# 6     0    0    0 DOWN    0   UP    0   UP    0    0 DOWN   0     2       6
# 7     0    0    0   UP   UP    0    0 DOWN    0   UP    0   0     4       3
# 8     0    0    0    0    0   UP   UP DOWN    0    0 DOWN   0     2       5
# 9     0    0   UP   UP    0    0    0    0    0   UP    0   0     7       0
# 10    0    0 DOWN    0    0    0    0 DOWN DOWN    0    0   0     2       4

print(paste("Up critical copper genes: ", nrow(r[which(r$numUp > r$numDown),]), sep = ""))
print(paste(r[which(r$numUp > r$numDown),"gene"], collapse=", "))
print(paste("Down critical copper genes: ", nrow(r[which(r$numUp < r$numDown),]), sep = ""))
print(paste(r[which(r$numUp < r$numDown),"gene"], collapse=", "))
print(paste("Other critical copper genes: ", nrow(r[which(r$numUp == r$numDown),]), sep = ""))
print(paste(r[which(r$numUp == r$numDown),"gene"], collapse=", "))
# > print(paste("Up critical copper genes: ", nrow(r[which(r$numUp > r$numDown),]), sep = ""))
# [1] "Up critical copper genes: 32"
# > print(paste(r[which(r$numUp > r$numDown),"gene"], collapse=", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, TP53, COA6, ARF1, GPC1, SORD, F5, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP, APC, ATOX1, ADAM17, STEAP3, MMGT1, AP1B1, CCND1, SPATA5, COMMD1, SUMF1, PARK7"
# > print(paste("Down critical copper genes: ", nrow(r[which(r$numUp < r$numDown),]), sep = ""))
# [1] "Down critical copper genes: 19"
# > print(paste(r[which(r$numUp < r$numDown),"gene"], collapse=", "))
# [1] "ALB, MAP1LC3A, SNCA, MAPT, JUN, CYP1A1, AOC3, PRNP, DBH, S100A12, AQP1, MT1X, XAF1, ANKRD9, LOXL1, AANAT, F8, HEPH, LCAT"
# > print(paste("Other critical copper genes: ", nrow(r[which(r$numUp == r$numDown),]), sep = ""))
# [1] "Other critical copper genes: 6"
# > print(paste(r[which(r$numUp == r$numDown),"gene"], collapse=", "))
# [1] "FOXO1, ANG, MT-CO2, BACE1, PRND, IL1A"

outDir <- paste(rootDir, "/Data/Output", sep = "")
write.csv(r,paste(outDir, '/criticalCopperGenesUpDown.csv', sep = ""), row.names = FALSE)

r <- r[which(r$frequency > 1),]
print(paste("Up critical copper genes: ", nrow(r[which(r$numUp > r$numDown),]), sep = ""))
print(paste(r[which(r$numUp > r$numDown),"gene"], collapse=", "))
print(paste("Down critical copper genes: ", nrow(r[which(r$numUp < r$numDown),]), sep = ""))
print(paste(r[which(r$numUp < r$numDown),"gene"], collapse=", "))
print(paste("Other critical copper genes: ", nrow(r[which(r$numUp == r$numDown),]), sep = ""))
print(paste(r[which(r$numUp == r$numDown),"gene"], collapse=", "))
# > print(paste("Up critical copper genes: ", nrow(r[which(r$numUp > r$numDown),]), sep = ""))
# [1] "Up critical copper genes: 21"
# > print(paste(r[which(r$numUp > r$numDown),"gene"], collapse=", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, TP53, COA6, ARF1, GPC1, SORD, F5, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"
# > print(paste("Down critical copper genes: ", nrow(r[which(r$numUp < r$numDown),]), sep = ""))
# [1] "Down critical copper genes: 13"
# > print(paste(r[which(r$numUp < r$numDown),"gene"], collapse=", "))
# [1] "ALB, MAP1LC3A, SNCA, MAPT, JUN, CYP1A1, AOC3, PRNP, DBH, S100A12, AQP1, MT1X, XAF1"
# > print(paste("Other critical copper genes: ", nrow(r[which(r$numUp == r$numDown),]), sep = ""))
# [1] "Other critical copper genes: 6"
# > print(paste(r[which(r$numUp == r$numDown),"gene"], collapse=", "))
# [1] "FOXO1, ANG, MT-CO2, BACE1, PRND, IL1A"

#================================================================

#================================================================
# (8) Obtain survival data, log2 data
# Run in server
#================================================================

# Obtain survival data, log data
obtainSurData=function(geneList) {
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- unique(subtypes$cancer.type)
  subtypes <- c(subtypes, "PAAD")
  subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
  subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
  n <- length(subtypes)
  surdata <- NA
  for (iType in 1:n) {
    cancertype <- subtypes[iType]
    outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
    
    # B. Prepare PAN.log2.counts dataframe
    
    # save(PAN.log2.counts.unique, file=paste(outDir, '/PAN.log2.counts.unique.Rdata', sep = ""))
    load(paste(outDir, '/PAN.log2.counts.unique.Rdata', sep = ""))  
    
    # create a subset of PAN.log2.counts that contains only the copper.critical
    PAN.log2.copper.critical <- PAN.log2.counts.unique[PAN.log2.counts.unique$gene %in% geneList, ]
    PAN.log2.copper.critical[1:5,1:3]
    nrow(PAN.log2.copper.critical)
    ncol(PAN.log2.copper.critical)
    # >     PAN.log2.copper.critical[1:5,1:3]
    # gene TCGA.OR.A5K3.01 TCGA.OR.A5J2.01
    # 32    AANAT          1.0000          0.0000
    # 4552 ADAM10          9.8281         12.6810
    # 4556 ADAM17          7.0799         10.9647
    # 5893    ALB          2.5850          4.3923
    # 6055    ANG          4.8851          8.4852
    # >     nrow(PAN.log2.copper.critical)
    # [1] 56
    # >     ncol(PAN.log2.copper.critical)
    # [1] 77
    PAN.log2.copper.critical <- PAN.log2.copper.critical[order(PAN.log2.copper.critical$gene, decreasing = FALSE),]
    PAN.log2.copper.critical[1:5,1:3]
    # >     PAN.log2.copper.critical[1:5,1:3]
    # gene TCGA.OR.A5K3.01 TCGA.OR.A5J2.01
    # 32    AANAT          1.0000          0.0000
    # 4552 ADAM10          9.8281         12.6810
    # 4556 ADAM17          7.0799         10.9647
    # 5893    ALB          2.5850          4.3923
    # 6055    ANG          4.8851          8.4852
    
    # C. Need to recreate the coxdata table
    
    if(cancertype != "PAAD") {
      x <- paste(rootDir, "/Data/Survival_SupplementalTable_S1_20171025_xena_sp", sep ="")
      PAN.surv <- read.table(x, header = TRUE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
      PAN.surv <- PAN.surv[, c("sample","OS", "X_PATIENT", "OS.time")]
      head(PAN.surv)
      nrow(PAN.surv)
      # >       head(PAN.surv)
      # sample OS    X_PATIENT OS.time
      # 1 TCGA-OR-A5J1-01  1 TCGA-OR-A5J1    1355
      # 2 TCGA-OR-A5J2-01  1 TCGA-OR-A5J2    1677
      # 3 TCGA-OR-A5J3-01  0 TCGA-OR-A5J3    2091
      # 4 TCGA-OR-A5J4-01  1 TCGA-OR-A5J4     423
      # 5 TCGA-OR-A5J5-01  1 TCGA-OR-A5J5     365
      # 6 TCGA-OR-A5J6-01  0 TCGA-OR-A5J6    2703
      # >       nrow(PAN.surv)
      # [1] 12591
    } else {
      x <- paste(rootDir, "/Data/TCGA-PAAD.survival.tsv", sep ="")
      PAN.surv <- read.table(x, header = TRUE, stringsAsFactors = FALSE)
      PAN.surv <- PAN.surv[, c("sample","OS", "X_PATIENT", "OS.time")]
      head(PAN.surv)
      nrow(PAN.surv)
      # > head(PAN.surv)
      # sample OS    X_PATIENT OS.time
      # 1 TCGA-HZ-7922-01A  0 TCGA-HZ-7922       4
      # 2 TCGA-HZ-A8P1-01A  0 TCGA-HZ-A8P1       7
      # 3 TCGA-IB-AAUM-01A  0 TCGA-IB-AAUM       8
      # 4 TCGA-RL-AAAS-01A  0 TCGA-RL-AAAS       9
      # 5 TCGA-US-A77G-11A  1 TCGA-US-A77G      12
      # 6 TCGA-US-A77G-01A  1 TCGA-US-A77G      12
      # > nrow(PAN.surv)
      # [1] 222
    }
    
    PAN.log2.copper.critical.T <- as.data.frame(t(PAN.log2.copper.critical)) #transpose data frame
    colnames(PAN.log2.copper.critical.T) <- PAN.log2.copper.critical$gene # assign gene symbol as colnames()
    PAN.log2.copper.critical.T <- PAN.log2.copper.critical.T[-1, ] # removes gene symbols as row values
    nGenes <- ncol(PAN.log2.copper.critical.T)
    PAN.log2.copper.critical.T[1:5,1:6]
    nrow(PAN.log2.copper.critical.T)
    ncol(PAN.log2.copper.critical.T)
    # >     PAN.log2.copper.critical.T[1:5,1:6]
    # AANAT  ADAM10  ADAM17     ALB     ANG  ANKRD9
    # TCGA.OR.A5K3.01  1.0000  9.8281  7.0799  2.5850  4.8851  9.4477
    # TCGA.OR.A5J2.01  0.0000 12.6810 10.9647  4.3923  8.4852 10.7712
    # TCGA.OR.A5LN.01  2.0000 10.0238  6.6852  6.9887  4.9383  8.5894
    # TCGA.OR.A5KY.01  3.5850 11.6443 10.1705  3.5850  8.8620 10.4208
    # TCGA.OR.A5LG.01  1.5850 11.8823  8.2943  2.5850  6.3045 10.4893
    # >     nrow(PAN.log2.copper.critical.T)
    # [1] 76
    # >     ncol(PAN.log2.copper.critical.T)
    # [1] 56
    
    if(nGenes != length(geneList)) {
      print("Please check this cancer type:")
      print(cancertype)
    }
    
    # convert counts characters into numeric
    i <- sapply(PAN.log2.copper.critical.T, is.character)
    PAN.log2.copper.critical.T[i] <- lapply(PAN.log2.copper.critical.T[i], as.numeric)
    
    # convert rownames in counts dataframe into a column
    PAN.log2.copper.critical.T$sample <- rownames(PAN.log2.copper.critical.T)
    PAN.log2.copper.critical.T <- PAN.log2.copper.critical.T %>% dplyr::select(sample, everything())
    
    # replace - with . to match the colnames in the counts dataframe
    tmp <- str_replace_all(PAN.surv$sample, "-", ".")
    PAN.surv$sample_ID <- tmp
    head(PAN.surv)
    # >     head(PAN.surv)
    # sample OS    X_PATIENT OS.time       sample_ID
    # 1 TCGA-OR-A5J1-01  1 TCGA-OR-A5J1    1355 TCGA.OR.A5J1.01
    # 2 TCGA-OR-A5J2-01  1 TCGA-OR-A5J2    1677 TCGA.OR.A5J2.01
    # 3 TCGA-OR-A5J3-01  0 TCGA-OR-A5J3    2091 TCGA.OR.A5J3.01
    # 4 TCGA-OR-A5J4-01  1 TCGA-OR-A5J4     423 TCGA.OR.A5J4.01
    # 5 TCGA-OR-A5J5-01  1 TCGA-OR-A5J5     365 TCGA.OR.A5J5.01
    # 6 TCGA-OR-A5J6-01  0 TCGA-OR-A5J6    2703 TCGA.OR.A5J6.01
    PAN.surv.tmp <- PAN.surv[ , -c(1,3)]
    head(PAN.surv.tmp)
    # >     head(PAN.surv.tmp)
    # OS OS.time       sample_ID
    # 1  1    1355 TCGA.OR.A5J1.01
    # 2  1    1677 TCGA.OR.A5J2.01
    # 3  0    2091 TCGA.OR.A5J3.01
    # 4  1     423 TCGA.OR.A5J4.01
    # 5  1     365 TCGA.OR.A5J5.01
    # 6  0    2703 TCGA.OR.A5J6.01
    
    PAN.log2.coxdata <- merge(PAN.log2.copper.critical.T, PAN.surv.tmp, by.x = "sample", by.y = "sample_ID")
    PAN.log2.coxdata <- PAN.log2.coxdata %>% dplyr::select((nGenes+2):(nGenes+3), everything())
    PAN.log2.coxdata <- PAN.log2.coxdata %>% dplyr::select(sample, everything())
    PAN.log2.coxdata[1:5,1:6]
    nrow(PAN.log2.coxdata)
    ncol(PAN.log2.coxdata)
    # > PAN.log2.coxdata[1:5,1:6]
    # sample OS OS.time  AANAT  ADAM10  ADAM17
    # 1 TCGA.OR.A5J1.01  1    1355 2.8074 11.8138  9.0570
    # 2 TCGA.OR.A5J2.01  1    1677 0.0000 12.6810 10.9647
    # 3 TCGA.OR.A5J3.01  0    2091 2.0000 10.9106  9.3005
    # 4 TCGA.OR.A5J5.01  1     365 2.8074  9.8872  8.2378
    # 5 TCGA.OR.A5J6.01  0    2703 2.0000 11.3756  9.3058
    # >     nrow(PAN.log2.coxdata)
    # [1] 76
    # >     ncol(PAN.log2.coxdata)
    # [1] 59
    
    #convert first column to rowname
    rownames(PAN.log2.coxdata) <- PAN.log2.coxdata[,1]
    PAN.log2.coxdata[,1] <- NULL
    PAN.log2.coxdata[1:5,1:6]
    # >     PAN.log2.coxdata[1:5,1:6]
    # OS OS.time  AANAT  ADAM10  ADAM17    ALB
    # TCGA.OR.A5J1.01  1    1355 2.8074 11.8138  9.0570 3.0000
    # TCGA.OR.A5J2.01  1    1677 0.0000 12.6810 10.9647 4.3923
    # TCGA.OR.A5J3.01  0    2091 2.0000 10.9106  9.3005 4.0000
    # TCGA.OR.A5J5.01  1     365 2.8074  9.8872  8.2378 8.6795
    # TCGA.OR.A5J6.01  0    2703 2.0000 11.3756  9.3058 3.1699
    
    # Univariate COX Analysis
    
    coxdata <- as.data.frame(PAN.log2.coxdata) # TCGA test data
    coxdata$cancertype <- cancertype
    coxdata <- coxdata %>% dplyr::select(cancertype, everything())
    
    if(iType == 1) {
      surdata <- coxdata
    } else {
      surdata <- rbind(surdata, coxdata)
    }
    
  }
  
  outDir <- paste(rootDir, "/Data/Output", sep = "")
  save(surdata, file=paste(outDir, '/surdata.Rdata', sep = ""))
  
  return(surdata)
}

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
geneList <- x$gene
surdata <- obtainSurData(geneList) # save in surdata.Rdata for 57 critical copper genes
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

table(surdata$cancertype)
# > table(surdata$cancertype)
# 
# ACC  AML BLCA BRCA COAD ESCA  GBM KICH KIRC KIRP  LGG LIHC LUAD LUSC OVCA PAAD
# 76  172  144 1211  277  179  165   91  508  189  348  239  274  192  418  153
# PCPG PRAD SKCM STAD THCA UCEC  UCS
# 179  373  331  409  560  182   57

#================================================================

#================================================================
# (9) Survival analysis
# Run in server
#================================================================

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 59

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# Read lasso genes
load(file = paste(outDir, "/lasso_min.Rdata", sep = "")) # return lasso_min
nrow(lasso_min)
head(lasso_min)
# > nrow(lasso_min)
# [1] 32
# > head(lasso_min)
# Active.Index Active.Coefficients lasso_gene
# 1            1          0.05136781     ADAM10
# 2            2          0.12289435        ANG
# 3            3         -0.02669716       AOC3
# 4            4         -0.01338193      AP1S1
# 5            5         -0.16595332        APP
# 6            6          0.02119514       AQP1

# All 56 critical copper genes
geneList56 <- x$gene
f56 <- paste(outDir, "/sur56.png", sep = "")
f56T <- paste(outDir, "/sur56T.png", sep = "")

# 39 popular critical cooper genes
geneList39 <- x[which(x$frequency > 1),]$gene
f39 <- paste(outDir, "/sur39.png", sep = "")
f39T <- paste(outDir, "/sur39T.png", sep = "")

# 31 up critical copper genes
geneList31 <- y[which(y$numUp > y$numDown),]$gene
f31 <- paste(outDir, "/sur31.png", sep = "")
f31T <- paste(outDir, "/sur31T.png", sep = "")

# 19 down critical copper genes
geneList19 <- y[which(y$numUp < y$numDown),]$gene
f19 <- paste(outDir, "/sur19.png", sep = "")
f19T <- paste(outDir, "/sur19T.png", sep = "")

# 20 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList20 <- y[which(y$numUp > y$numDown),]$gene
f20 <- paste(outDir, "/sur20.pdf", sep = "")

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene
f13 <- paste(outDir, "/sur13.pdf", sep = "")

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# > head(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# > nrow(uni_gene)
# [1] 35

# 32 lasso genes
# 18 up lasso genes
geneList18 <- geneList20[which(geneList20 %in% lasso_min$lasso_gene)]
f18 <- paste(outDir, "/sur18.pdf", sep = "")
length(geneList18)
print(paste(geneList18, collapse = ", "))
# > length(geneList18)
# [1] 18
# > print(paste(geneList18, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, ARF1, GPC1, SORD, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"

# 10 down lasso genes
geneList10 <- geneList13[which(geneList13 %in% lasso_min$lasso_gene)]
f10 <- paste(outDir, "/sur10.pdf", sep = "")
length(geneList10)
print(paste(geneList10, collapse = ", "))
# > length(geneList10)
# [1] 10
# > print(paste(geneList10, collapse = ", "))
# [1] "MAP1LC3A, SNCA, MAPT, CYP1A1, AOC3, PRNP, S100A12, AQP1, MT1X, XAF1"

# Cox genes
geneList18 <- geneList20[which(geneList20 %in% uni_gene$x)]
f18 <- paste(outDir, "/sur18.pdf", sep = "")
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]
f12 <- paste(outDir, "/sur12.pdf", sep = "")
surAnalysis(surdata, geneList18, numGroup=2, K=20, alpha=0.5, f18, w = 9, h = 6)
surAnalysis(surdata, geneList12, numGroup=2, K=20, alpha=0.5, f12, w = 9, h = 6)
print(paste(geneList18, collapse = ", "))
print(paste(geneList12, collapse = ", "))
# > print(paste(geneList18, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, ARF1, GPC1, SORD, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"
# > print(paste(geneList12, collapse = ", "))
# [1] "MAP1LC3A, SNCA, MAPT, JUN, CYP1A1, AOC3, PRNP, DBH, S100A12, AQP1, MT1X, XAF1"

surAnalysis(surdata, geneList56, numGroup=2, K=20, alpha=0.5, f56, w = 800, h = 600)
surAnalysis(surdata, geneList39, numGroup=2, K=20, alpha=0.5, f39, w = 800, h = 600)
surAnalysis(surdata, geneList31, numGroup=2, K=20, alpha=0.5, f31, w = 800, h = 600)
surAnalysis(surdata, geneList19, numGroup=2, K=20, alpha=0.5, f19, w = 800, h = 600)
surAnalysis(surdata, geneList20, numGroup=2, K=20, alpha=0.5, f20, w = 8, h = 6)
surAnalysis(surdata, geneList13, numGroup=2, K=20, alpha=0.5, f13, w = 8, h = 6)
surAnalysis(surdata, geneList18, numGroup=2, K=20, alpha=0.5, f18, w = 8.5, h = 6)
surAnalysis(surdata, geneList10, numGroup=2, K=20, alpha=0.5, f10, w = 8.5, h = 6)

# # Survival analysis
# surAnalysis(surdata, geneList56, numGroup=2, K=20, alpha=0.5, f56, w = 800, h = 600)
# surAnalysis(surdata, geneList56, numGroup=3, K=20, alpha=0.5, f56T, w = 800, h = 600)
# surAnalysis(surdata, geneList39, numGroup=2, K=20, alpha=0.5, f39, w = 800, h = 600)
# surAnalysis(surdata, geneList39, numGroup=3, K=20, alpha=0.5, f39T, w = 800, h = 600)
# surAnalysis(surdata, geneList31, numGroup=2, K=20, alpha=0.5, f31, w = 800, h = 600)
# surAnalysis(surdata, geneList31, numGroup=3, K=20, alpha=0.5, f31T, w = 800, h = 600)
# surAnalysis(surdata, geneList19, numGroup=2, K=20, alpha=0.5, f19, w = 800, h = 600)
# surAnalysis(surdata, geneList19, numGroup=3, K=20, alpha=0.5, f19T, w = 800, h = 600)

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
n <- length(subtypes)

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
criticalCopperGenes <- read.csv(fileName)
criticalCopperGenes[1:5,1:6]
ncol(criticalCopperGenes)
colnames(criticalCopperGenes)
# > criticalCopperGenes[1:5,1:6]
# gene frequency ACC AML BLCA BRCA
# 1     CDK1        18   1   1    1    1
# 2      ALB        13   0   1    0    1
# 3    AP1S1        11   0   1    0    0
# 4    CASP3         9   0   0    0    1
# 5 MAP1LC3A         8   0   1    0    1
# > ncol(criticalCopperGenes)
# [1] 26
# > colnames(criticalCopperGenes)
# [1] "gene"      "frequency" "ACC"       "AML"       "BLCA"      "BRCA"
# [7] "COAD"      "ESCA"      "GBM"       "KICH"      "KIRC"      "KIRP"
# [13] "LGG"       "LIHC"      "LUAD"      "LUSC"      "OVCA"      "PAAD"
# [19] "PCPG"      "PRAD"      "SKCM"      "STAD"      "THCA"      "UCEC"
# [25] "UCS"       "summary"

for (iType in 1:n) {
  cancertype <- subtypes[iType]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  dat <- surdata[which(surdata$cancertype == cancertype),]
  
  # Critical copper genes
  geneList <- criticalCopperGenes[which(criticalCopperGenes[,iType+2] == 1),]$gene
  f <- paste(outDir, "/", cancertype, "Sur.pdf", sep = "")
  fT <- paste(outDir, "/", cancertype, "SurT.pdf", sep = "")
  
  surAnalysis(dat, geneList, numGroup=2, K=20, alpha=0.5, f, w = 9, h = 6)
  surAnalysis(dat, geneList, numGroup=3, K=20, alpha=0.5, fT, w = 9, h = 6)
}

#================================================================

#================================================================
# (10) Survival analysis
# Univariate Cox Analysis (Katana)
# Run in server
#================================================================
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
criticalCopperGenes <- read.csv(fileName)
n <- nrow(criticalCopperGenes)
criticalCopperGenes[1:5,1:6]
n
# > criticalCopperGenes[1:5,1:6]
# gene frequency ACC AML BLCA BRCA
# 1     CDK1        18   1   1    1    1
# 2      ALB        13   0   1    0    1
# 3    AP1S1        11   0   1    0    0
# 4    CASP3         9   0   0    0    1
# 5 MAP1LC3A         8   0   1    0    1
# > n
# [1] 57

r <- criticalCopperGenes[,-ncol(criticalCopperGenes)] # remove summary
r[,3:(ncol(criticalCopperGenes)-1)] <- 2

for (i in 3:(ncol(criticalCopperGenes)-1)) {
  cancertype <- colnames(criticalCopperGenes)[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(paste(outDir, '/unicox.Rdata', sep = ""))
  idx <- nrow(unicox)
  for (j in 1:idx) {
    gene <- unicox$Variable[j]
    r[which(r$gene == gene),i] <- unicox$global.pval[j]
  }
}

# check r
head(r)
# > head(r)
# gene frequency      ACC    AML   BLCA    BRCA   COAD  ESCA    GBM
# 1     CDK1        18 3.87e-09 0.6490 0.0569 0.01180 0.2390 0.623 2.0000
# 2      ALB        13 2.00e+00 0.3390 2.0000 0.01580 0.2610 2.000 0.3970
# 3    AP1S1        11 2.00e+00 0.3730 2.0000 2.00000 0.8110 2.000 2.0000
# 4    CASP3         9 2.00e+00 2.0000 2.0000 0.76000 0.2750 2.000 0.0849
# 5 MAP1LC3A         8 2.00e+00 0.0404 2.0000 0.07150 2.0000 0.598 2.0000
# 6     SNCA         8 2.00e+00 2.0000 2.0000 0.00972 0.0326 2.000 0.1290
# KICH KIRC   KIRP      LGG  LIHC  LUAD   LUSC   OVCA   PAAD     PCPG  PRAD
# 1 3.42e-05    2 2.0000 9.33e-07 0.114 0.111 2.0000 0.1640 0.0104 0.004730 0.272
# 2 8.76e-01    2 0.5070 6.30e-04 2.000 2.000 2.0000 0.2880 0.1450 2.000000 0.334
# 3 2.00e+00    2 0.0469 2.00e+00 0.550 0.958 2.0000 0.7110 2.0000 0.846000 0.203
# 4 9.26e-01    2 0.3560 4.91e-06 2.000 2.000 2.0000 0.0423 0.7510 2.000000 2.000
# 5 2.00e+00    2 2.0000 1.21e-01 2.000 2.000 0.0298 2.0000 2.0000 0.000304 2.000
# 6 2.00e+00    2 2.0000 5.76e-01 2.000 2.000 2.0000 0.1550 2.0000 0.934000 2.000
# SKCM   STAD  THCA  UCEC   UCS
# 1 0.836 0.4780 0.898 0.311 2.000
# 2 0.721 0.0395 0.049 2.000 2.000
# 3 2.000 0.0245 2.000 0.494 0.572
# 4 2.000 2.0000 2.000 0.526 2.000
# 5 0.634 2.0000 0.147 2.000 2.000
# 6 0.674 2.0000 2.000 0.703 2.000

rr <- r
rr$Total <- rowSums(r[,-c(1,2)] < 2)
rr$NumSig <- rowSums(r[,-c(1,2)] <= 0.05)
rr <- rr[order(rr$NumSig, decreasing = TRUE),]
rr$gene[1:10]
# > rr$gene[1:10]
# [1] "CDK1"     "ALB"      "MAP1LC3A" "TMPRSS6"  "AP1S1"    "CASP3"
# [7] "SNCA"     "CYP1A1"   "XIAP"     "AOC3"

outDir <- paste(rootDir, "/Data/Output", sep = "")
write.csv(rr, paste(outDir, '/pCox.csv', sep = ""), row.names = FALSE)

# Draw the chart
dat <- read.csv(paste(outDir, '/pCox.csv', sep = ""))
#dat <- dat[1:10,]
# Get only genes which are significant in at least one cancer type
dat <- dat[which(dat$NumSig > 0),]

rownames(dat) <- dat$gene
n <- ncol(dat)
dat <- dat[,-c(1,2,n-1,n)]
dat <- as.matrix(dat)

tmp <- dat
tmp <- ifelse(tmp == 2, NA, tmp)
tmp2 <- melt(tmp)
colnames(tmp2) <- c("Gene", "Cancer", "pvalue")
tmp2$Label <- ifelse(tmp2$pvalue > 0.05, "", "*")

# tmp2 %>%
#   ggplot(aes(Gene, Cancer, fill = pvalue)) + geom_tile() +
#   geom_text(aes(label = Label)) + 
#   theme(axis.title = element_blank())

# tmp3 <- tmp2 %>%
#   arrange(pvalue) %>%
#   mutate(y = fct_reorder(Gene, pvalue, count, .desc = FALSE),
#          y1 = as.integer(Gene),
#          x = factor(Cancer),
#          x1 = as.integer(Cancer))

tmpx <- tmp2
nRow <- nrow(tmpx)
tmpx$orderCol <- 0
for (i in 1:nRow) {
  gene <- tmpx$Gene[i]
  tmpxx <- tmpx[which(tmpx$Gene == gene),]
  tmpxx <- tmpxx[which(!is.na(tmpxx$pvalue)),]
  tmpxx <- tmpxx[which(tmpxx$pvalue <= 0.05),]
  tmpx$orderCol[i] <- nrow(tmpxx)  
}
tmp2 <- tmpx

# Prepare data
tmp3 <- tmp2
tmp3 <- arrange(tmp3, orderCol)
tmp3 <- mutate(tmp3, y = fct_reorder(Gene, orderCol, min, .desc = FALSE),
               y1 = as.integer(y),
               x = factor(Cancer),
               x1 = as.integer(x))

labels_y <- levels(tmp3$y)
breaks_y <- seq_along(labels_y)

labels_x <- levels(tmp3$x)
breaks_x <- seq_along(labels_x)

outFile <- paste(outDir, "/uniCox.pdf", sep = "")
pdf(file = outFile,  width = 9, height = 6, onefile=FALSE)
Note <- c("pvalue < 0.05")
ggplot(tmp3, aes(x=x1, y=y1))+
  geom_tile(aes(fill=pvalue), color="white", size=0.25) + 
  scale_y_continuous(breaks = breaks_y, labels = labels_y) +
  scale_x_continuous(breaks = breaks_x, labels = labels_x) +
  scale_fill_continuous(low="#0066E7", high="#e5effc", 
                        guide="colorbar", na.value="#CCCCCC") +
  theme(panel.grid.major = element_blank(), axis.title=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=15)) +
  geom_text(aes(label = Label, colour=Note)) +
  scale_colour_manual(values=c("red"))
dev.off()
#================================================================

#================================================================
# (11) Survival analysis
# Univariate Cox Analysis (Katana) for all cancer types
# Run in server
#================================================================

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
x[1:5,1:4]
nrow(x)
ncol(x)
y[1:5,1:4]
nrow(y)
ncol(y)
colnames(y)
# > x[1:5,1:4]
# gene frequency ACC AML
# 1     CDK1        18   1   1
# 2      ALB        13   0   1
# 3    AP1S1        11   0   1
# 4    CASP3         9   0   0
# 5 MAP1LC3A         8   0   1
# > nrow(x)
# [1] 57
# > ncol(x)
# [1] 26
# > y[1:5,1:4]
# gene frequency ACC  AML
# 1     CDK1        18  UP   UP
# 2      ALB        13   0 DOWN
# 3    AP1S1        11   0 DOWN
# 4    CASP3         9   0    0
# 5 MAP1LC3A         8   0 DOWN
# > nrow(y)
# [1] 57
# > ncol(y)
# [1] 27
# > colnames(y)
# [1] "gene"      "frequency" "ACC"       "AML"       "BLCA"      "BRCA"
# [7] "COAD"      "ESCA"      "GBM"       "KICH"      "KIRC"      "KIRP"
# [13] "LGG"       "LIHC"      "LUAD"      "LUSC"      "OVCA"      "PAAD"
# [19] "PCPG"      "PRAD"      "SKCM"      "STAD"      "THCA"      "UCEC"
# [25] "UCS"       "numUp"     "numDown"

# All 57 critical copper genes
geneList57 <- x$gene
f57 <- paste(outDir, "/uni57.csv", sep = "")

# 32 up critical copper genes
geneList32 <- y[which(y$numUp > y$numDown),]$gene
f32 <- paste(outDir, "/uni32.csv", sep = "")

# 19 down critical copper genes
geneList19 <- y[which(y$numUp < y$numDown),]$gene
f19 <- paste(outDir, "/uni19.csv", sep = "")

# 40 popular critical cooper genes
geneList40 <- x[which(x$frequency > 1),]$gene
f40 <- paste(outDir, "/uni40.csv", sep = "")

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene
f21 <- paste(outDir, "/uni21.csv", sep = "")

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene
f13 <- paste(outDir, "/uni13.csv", sep = "")

# Univariate Cox Analysis must use log2(count+1) data

# A. First select copper genes from critical list

# geneList40
# geneList21
# geneList13

# B. Prepare PAN.log2.counts dataframe

# surdata

# C. Need to recreate the coxdata table

# surdata

####----First Univariate COX analysis to identify genes----####

# Univariate COX Analysis

coxdata <- as.data.frame(surdata)
coxdata[1:5,1:4]
nrow(coxdata)
ncol(coxdata)
# > coxdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(coxdata)
# [1] 6727
# > ncol(coxdata)
# [1] 60

# to check distibution of gene expression (Optional) - they should have a normal distribution
# check gene MT1X
#ggplot(coxdata, aes(x = MT1X)) + geom_histogram(color = "black", fill = "white") 

res_mod = ezcox(coxdata, time = "OS.time", status = "OS", covariates = geneList40, global_method = c("likelihood", "wald", "logrank"), return_models = TRUE) # for TCGA test data
mds = get_models(res_mod)
str(mds, max.level = 1)
show_models(mds)

cox_res <- ezcox(coxdata, time = "OS.time", status = "OS", covariates = geneList40, global_method = c("likelihood", "wald", "logrank"))
unicox <- cox_res
unicox[1:5,1:6]
nrow(unicox)
ncol(unicox)
# > unicox[1:5,1:6]
# # A tibble: 5 × 6
# Variable is_control contrast_level ref_level n_contrast n_ref
# <chr>    <lgl>      <chr>          <chr>          <int> <int>
#   1 CDK1     FALSE      CDK1           CDK1            6678  6678
# 2 ALB      FALSE      ALB            ALB             6678  6678
# 3 AP1S1    FALSE      AP1S1          AP1S1           6678  6678
# 4 CASP3    FALSE      CASP3          CASP3           6678  6678
# 5 MAP1LC3A FALSE      MAP1LC3A       MAP1LC3A        6678  6678
# > nrow(unicox)
# [1] 40
# > ncol(unicox)
# [1] 12

write.csv(unicox,file = paste(outDir, "/unicox.csv", sep =""),row.names = FALSE)
save(unicox, file=paste(outDir, '/unicox.Rdata', sep = ""))

save(coxdata, file=paste(outDir, '/coxdata.Rdata', sep = ""))

save(res_mod, file=paste(outDir, '/res_mod.Rdata', sep = ""))

######Get significant genes#######

load(paste(outDir, '/unicox.Rdata', sep = ""))
unicox[1:5,1:6]
nrow(unicox)
ncol(unicox)
# > unicox[1:5,1:6]
# # A tibble: 5 × 6
# Variable is_control contrast_level ref_level n_contrast n_ref
# <chr>    <lgl>      <chr>          <chr>          <int> <int>
#   1 CDK1     FALSE      CDK1           CDK1            6678  6678
# 2 ALB      FALSE      ALB            ALB             6678  6678
# 3 AP1S1    FALSE      AP1S1          AP1S1           6678  6678
# 4 CASP3    FALSE      CASP3          CASP3           6678  6678
# 5 MAP1LC3A FALSE      MAP1LC3A       MAP1LC3A        6678  6678
# > nrow(unicox)
# [1] 40
# > ncol(unicox)
# [1] 12

unicox_dif <- unicox[which(unicox$global.pval < 0.05), ]
nrow(unicox_dif)
# > nrow(unicox_dif)
# [1] 35
uni_gene <- as.character(unicox_dif$Variable)
uni_gene
# > uni_gene
# [1] "CDK1"     "AP1S1"    "CASP3"    "MAP1LC3A" "SNCA"     "TMPRSS6"
# [7] "MAPT"     "GSK3B"    "JUN"      "APP"      "CYP1A1"   "COX17"
# [13] "XIAP"     "FOXO1"    "ARF1"     "GPC1"     "AOC3"     "SORD"
# [19] "PRNP"     "ATP7A"    "SP1"      "MT-CO1"   "DBH"      "SLC11A2"
# [25] "ANG"      "S100A12"  "ATP6AP1"  "ADAM10"   "MT-CO2"   "CP"
# [31] "PRND"     "AQP1"     "MT1X"     "IL1A"     "XAF1"

print(paste(uni_gene, collapse=", "))
# > print(paste(uni_gene, collapse=", "))
# [1] "CDK1, AP1S1, CASP3, MAP1LC3A, SNCA, TMPRSS6, MAPT, GSK3B, JUN, APP, CYP1A1, COX17, XIAP, FOXO1, ARF1, GPC1, AOC3, SORD, PRNP, ATP7A, SP1, MT-CO1, DBH, SLC11A2, ANG, S100A12, ATP6AP1, ADAM10, MT-CO2, CP, PRND, AQP1, MT1X, IL1A, XAF1"

write.csv(uni_gene,file = paste(outDir, "/prognostic_genes.csv", sep =""),row.names = FALSE)

#======================
######constructing risk-score system, LASSO model#######

load(paste(outDir, '/coxdata.Rdata', sep = ""))
coxdata[1:5,1:6]
nrow(coxdata)
ncol(coxdata)
colnames(coxdata)
# > coxdata[1:5,1:6]
# cancertype OS OS.time  AANAT  ADAM10  ADAM17
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074 11.8138  9.0570
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000 12.6810 10.9647
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000 10.9106  9.3005
# TCGA.OR.A5J5.01        ACC  1     365 2.8074  9.8872  8.2378
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000 11.3756  9.3058
# > nrow(coxdata)
# [1] 6727
# > ncol(coxdata)
# [1] 59
# > colnames(coxdata)
# [1] "cancertype" "OS"         "OS.time"    "AANAT"      "ADAM10"
# [6] "ADAM17"     "ALB"        "ANG"        "ANKRD9"     "AOC3"
# [11] "AP1B1"      "AP1S1"      "APC"        "APP"        "AQP1"
# [16] "ARF1"       "ATOX1"      "ATP6AP1"    "ATP7A"      "BACE1"
# [21] "CASP3"      "CCND1"      "CDK1"       "COMMD1"     "COX17"
# [26] "CP"         "CYP1A1"     "DBH"        "F5"         "F8"
# [31] "FOXO1"      "GPC1"       "GSK3B"      "HEPH"       "IL1A"
# [36] "JUN"        "LCAT"       "LOXL1"      "MAP1LC3A"   "MAPT"
# [41] "MMGT1"      "MT-CO1"     "MT-CO2"     "MT1X"       "PARK7"
# [46] "PRND"       "PRNP"       "S100A12"    "SLC11A2"    "SNCA"
# [51] "SORD"       "SP1"        "SPATA5"     "STEAP3"     "SUMF1"
# [56] "TMPRSS6"    "TP53"       "XAF1"       "XIAP"

ncol <- ncol(coxdata)
irg_expr <- coxdata[ ,c(4:ncol)]

x <- as.matrix(irg_expr)
x <- x[,which(colnames(x) %in% unicox_dif$Variable)]
y <- data.matrix(Surv(coxdata$OS.time,coxdata$OS))
x[1:5, 1:4]
nrow(x)
head(y)
nrow(y)
# > x[1:5, 1:4]
# ADAM10    ANG    AOC3   AP1S1
# TCGA.OR.A5J1.01 11.8138 6.9844  8.0279 10.6627
# TCGA.OR.A5J2.01 12.6810 8.4852 10.3151 12.6517
# TCGA.OR.A5J3.01 10.9106 9.1847  9.2897 13.5362
# TCGA.OR.A5J5.01  9.8872 7.9170  9.0461 11.4528
# TCGA.OR.A5J6.01 11.3756 8.1944 12.6156 11.5196
# > nrow(x)
# [1] 6727
# > head(y)
# time status
# [1,] 1355      1
# [2,] 1677      1
# [3,] 2091      0
# [4,]  365      1
# [5,] 2703      0
# [6,]  490      1
# > nrow(y)
# [1] 6727

# remove time with NA or <= 0
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
# [1]   91   94  105  110  127  145  152  154  155  165  168  173  176  179  183
# [16]  185  186  187  195  196  198  199  200  219  380  567  595  597  606  617
# [31]  621  623  631  642 1150 1157 1167 1565 1580 1659 1660 1678 1729 1731 2284
# [46] 2285 2979 3018 3261 3262 3394 3395 3586 3604 3610 3801 3802 3805 3806 3807
# [61] 3808 3813 3814 3841 4004 4023 4024 4076 4230 5219 5231 5232 5247 5301 5326
# [76] 5330 5350 5435 5446 5461 5521 5538 5539 5540 5541 5542 5543 5544 5545 5546
# [91] 5547 5548 5552 5553 5587 5588 5631 5641 5684 5692 5696 5700 5702 5781 5828
# [106] 6414 6509 6544 6727
# > x[1:5, 1:4]
# ADAM10    ANG    AOC3   AP1S1
# TCGA.OR.A5J1.01 11.8138 6.9844  8.0279 10.6627
# TCGA.OR.A5J2.01 12.6810 8.4852 10.3151 12.6517
# TCGA.OR.A5J3.01 10.9106 9.1847  9.2897 13.5362
# TCGA.OR.A5J5.01  9.8872 7.9170  9.0461 11.4528
# TCGA.OR.A5J6.01 11.3756 8.1944 12.6156 11.5196
# > nrow(x)
# [1] 6618
# > head(y)
# time status
# [1,] 1355      1
# [2,] 1677      1
# [3,] 2091      0
# [4,]  365      1
# [5,] 2703      0
# [6,]  490      1
# > nrow(y)
# [1] 6618

# Coefficient profiles in the LASSO regression model.
pdf(file = paste(outDir, '/coefPro.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
fit0 <- glmnet(x, y, family = "cox", alpha = 1, nlambda = 1000)
plot(fit0)
dev.off()
pdf(file = paste(outDir, '/coefProLam.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
plot(fit0, xvar="lambda", label=FALSE)
dev.off()

set.seed(1)

# Cross-validation for tuning parameter screening in the LASSO regression model. 
pdf(file = paste(outDir, '/crossVal.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
cv.fit <- cv.glmnet(x, y,
                    family="cox",
                    maxit = 1000000,
                    alpha=1)
print(cv.fit)
plot(cv.fit)
dev.off()

# LASSO_gene table
pdf(file = paste(outDir, '/lasso_fit.pdf', sep = ""),  width = 5.5, height = 4, onefile=FALSE)
fit <- glmnet(x, y, alpha = 1, family='cox',lambda=cv.fit$lambda.min)
plot(fit)
coef(fit)
dev.off()
# > coef(fit)
# 35 x 1 sparse Matrix of class "dgCMatrix"
# s0
# ADAM10    0.051367811
# ANG       0.122894349
# AOC3     -0.026697160
# AP1S1    -0.013381933
# APP      -0.165953320
# AQP1      0.021195137
# ARF1     -0.196391961
# ATP6AP1  -0.018777436
# ATP7A    -0.053776030
# CASP3     0.047855806
# CDK1      0.222878440
# COX17    -0.352159897
# CP        0.040369618
# CYP1A1    0.030710229
# DBH       .
# FOXO1     0.045857538
# GPC1      0.142751617
# GSK3B     0.137378822
# IL1A      0.031853610
# JUN       .
# MAP1LC3A  0.089213469
# MAPT     -0.019495688
# MT-CO1    0.139414786
# MT-CO2    .
# MT1X     -0.025790189
# PRND     -0.021438062
# PRNP      0.025102772
# S100A12   0.060583996
# SLC11A2  -0.065724355
# SNCA      0.004206572
# SORD     -0.142404992
# SP1      -0.153113728
# TMPRSS6  -0.004168475
# XAF1     -0.001730944
# XIAP     -0.086420616

Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Coefficients
# > Coefficients
# 35 x 1 sparse Matrix of class "dgCMatrix"
# 1
# ADAM10    0.051367811
# ANG       0.122894349
# AOC3     -0.026697160
# AP1S1    -0.013381933
# APP      -0.165953320
# AQP1      0.021195137
# ARF1     -0.196391961
# ATP6AP1  -0.018777436
# ATP7A    -0.053776030
# CASP3     0.047855806
# CDK1      0.222878440
# COX17    -0.352159897
# CP        0.040369618
# CYP1A1    0.030710229
# DBH       .
# FOXO1     0.045857538
# GPC1      0.142751617
# GSK3B     0.137378822
# IL1A      0.031853610
# JUN       .
# MAP1LC3A  0.089213469
# MAPT     -0.019495688
# MT-CO1    0.139414786
# MT-CO2    .
# MT1X     -0.025790189
# PRND     -0.021438062
# PRNP      0.025102772
# S100A12   0.060583996
# SLC11A2  -0.065724355
# SNCA      0.004206572
# SORD     -0.142404992
# SP1      -0.153113728
# TMPRSS6  -0.004168475
# XAF1     -0.001730944
# XIAP     -0.086420616

Active.Index
Active.Coefficients
# > Active.Index
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 16 17 18 19 21 22 23 25 26 27 28
# [26] 29 30 31 32 33 34 35
# > Active.Coefficients
# [1]  0.051367811  0.122894349 -0.026697160 -0.013381933 -0.165953320
# [6]  0.021195137 -0.196391961 -0.018777436 -0.053776030  0.047855806
# [11]  0.222878440 -0.352159897  0.040369618  0.030710229  0.045857538
# [16]  0.142751617  0.137378822  0.031853610  0.089213469 -0.019495688
# [21]  0.139414786 -0.025790189 -0.021438062  0.025102772  0.060583996
# [26] -0.065724355  0.004206572 -0.142404992 -0.153113728 -0.004168475
# [31] -0.001730944 -0.086420616

lasso_gene <- row.names(Coefficients)[Active.Index]
lasso_min <- data.frame(Active.Index,Active.Coefficients,lasso_gene)
lasso_min
# > lasso_min
# Active.Index Active.Coefficients lasso_gene
# 1             1         0.051367811     ADAM10
# 2             2         0.122894349        ANG
# 3             3        -0.026697160       AOC3
# 4             4        -0.013381933      AP1S1
# 5             5        -0.165953320        APP
# 6             6         0.021195137       AQP1
# 7             7        -0.196391961       ARF1
# 8             8        -0.018777436    ATP6AP1
# 9             9        -0.053776030      ATP7A
# 10           10         0.047855806      CASP3
# 11           11         0.222878440       CDK1
# 12           12        -0.352159897      COX17
# 13           13         0.040369618         CP
# 14           14         0.030710229     CYP1A1
# 15           16         0.045857538      FOXO1
# 16           17         0.142751617       GPC1
# 17           18         0.137378822      GSK3B
# 18           19         0.031853610       IL1A
# 19           21         0.089213469   MAP1LC3A
# 20           22        -0.019495688       MAPT
# 21           23         0.139414786     MT-CO1
# 22           25        -0.025790189       MT1X
# 23           26        -0.021438062       PRND
# 24           27         0.025102772       PRNP
# 25           28         0.060583996    S100A12
# 26           29        -0.065724355    SLC11A2
# 27           30         0.004206572       SNCA
# 28           31        -0.142404992       SORD
# 29           32        -0.153113728        SP1
# 30           33        -0.004168475    TMPRSS6
# 31           34        -0.001730944       XAF1
# 32           35        -0.086420616       XIAP

print(paste(lasso_gene, collapse=", "))
# > print(paste(lasso_gene, collapse=", "))
# [1] "ADAM10, ANG, AOC3, AP1S1, APP, AQP1, ARF1, ATP6AP1, ATP7A, CASP3, CDK1, COX17, CP, CYP1A1, FOXO1, GPC1, GSK3B, IL1A, MAP1LC3A, MAPT, MT-CO1, MT1X, PRND, PRNP, S100A12, SLC11A2, SNCA, SORD, SP1, TMPRSS6, XAF1, XIAP"

save(lasso_min,file = paste(outDir, "/lasso_min.Rdata", sep = ""))
save(cv.fit,fit,lasso_gene,file = paste(outDir, "/Lasso_model_min.Rdata", sep = ""))
#=====================

#================================================================

#================================================================
# (12) Latex table for LASSO, pan cancer and LGG
#================================================================
# For pan cancer
outDir <- paste(rootDir, "/Data/Output", sep = "")
load(file = paste(outDir, "/lasso_min.Rdata", sep = "")) # return lasso_min
head(lasso_min)
nrow(lasso_min)
# > head(lasso_min)
# Active.Index Active.Coefficients lasso_gene
# 1            1          0.05136781     ADAM10
# 2            2          0.12289435        ANG
# 3            3         -0.02669716       AOC3
# 4            4         -0.01338193      AP1S1
# 5            5         -0.16595332        APP
# 6            6          0.02119514       AQP1
# > nrow(lasso_min)
# [1] 32

t <- lasso_min
t <- t[,c(3,2)]
colnames(t) <- c("Prognostic gene", "Coefficient")
head(t)
# > head(t)
# Prognostic gene Coefficient
# 1          ADAM10  0.05136781
# 2             ANG  0.12289435
# 3            AOC3 -0.02669716
# 4           AP1S1 -0.01338193
# 5             APP -0.16595332
# 6            AQP1  0.02119514

t[,2]<-format(as.numeric(t[,2]),big.mark=",")
print(latex.table.by(t, format.args = list(digits = 2, format = c("s","d"))), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)
# > print(latex.table.by(t, format.args = list(digits = 2, format = c("s","d"))), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)
# % latex table generated in R 4.0.2 by xtable 1.8-4 package
# % Tue Oct 11 12:03:03 2022
# \begin{table}[ht]
# \centering
# \begin{tabular}{|c|r|}
# \hline
# Prognostic gene & Coefficient \\
# \hline
# \cline{1-2}\multirow{ 1 }{*}{ ADAM10 } &  0.051367811 \\
# \cline{1-2}\multirow{ 1 }{*}{ ANG } &  0.122894349 \\
# \cline{1-2}\multirow{ 1 }{*}{ AOC3 } & -0.026697160 \\
# \cline{1-2}\multirow{ 1 }{*}{ AP1S1 } & -0.013381933 \\
# \cline{1-2}\multirow{ 1 }{*}{ APP } & -0.165953320 \\
# \cline{1-2}\multirow{ 1 }{*}{ AQP1 } &  0.021195137 \\
# \cline{1-2}\multirow{ 1 }{*}{ ARF1 } & -0.196391961 \\
# \cline{1-2}\multirow{ 1 }{*}{ ATP6AP1 } & -0.018777436 \\
# \cline{1-2}\multirow{ 1 }{*}{ ATP7A } & -0.053776030 \\
# \cline{1-2}\multirow{ 1 }{*}{ CASP3 } &  0.047855806 \\
# \cline{1-2}\multirow{ 1 }{*}{ CDK1 } &  0.222878440 \\
# \cline{1-2}\multirow{ 1 }{*}{ COX17 } & -0.352159897 \\
# \cline{1-2}\multirow{ 1 }{*}{ CP } &  0.040369618 \\
# \cline{1-2}\multirow{ 1 }{*}{ CYP1A1 } &  0.030710229 \\
# \cline{1-2}\multirow{ 1 }{*}{ FOXO1 } &  0.045857538 \\
# \cline{1-2}\multirow{ 1 }{*}{ GPC1 } &  0.142751617 \\
# \cline{1-2}\multirow{ 1 }{*}{ GSK3B } &  0.137378822 \\
# \cline{1-2}\multirow{ 1 }{*}{ IL1A } &  0.031853610 \\
# \cline{1-2}\multirow{ 1 }{*}{ MAP1LC3A } &  0.089213469 \\
# \cline{1-2}\multirow{ 1 }{*}{ MAPT } & -0.019495688 \\
# \cline{1-2}\multirow{ 1 }{*}{ MT-CO1 } &  0.139414786 \\
# \cline{1-2}\multirow{ 1 }{*}{ MT1X } & -0.025790189 \\
# \cline{1-2}\multirow{ 1 }{*}{ PRND } & -0.021438062 \\
# \cline{1-2}\multirow{ 1 }{*}{ PRNP } &  0.025102772 \\
# \cline{1-2}\multirow{ 1 }{*}{ S100A12 } &  0.060583996 \\
# \cline{1-2}\multirow{ 1 }{*}{ SLC11A2 } & -0.065724355 \\
# \cline{1-2}\multirow{ 1 }{*}{ SNCA } &  0.004206572 \\
# \cline{1-2}\multirow{ 1 }{*}{ SORD } & -0.142404992 \\
# \cline{1-2}\multirow{ 1 }{*}{ SP1 } & -0.153113728 \\
# \cline{1-2}\multirow{ 1 }{*}{ TMPRSS6 } & -0.004168475 \\
# \cline{1-2}\multirow{ 1 }{*}{ XAF1 } & -0.001730944 \\
# \cline{1-2}\multirow{ 1 }{*}{ XIAP } & -0.086420616 \\
# \hline
# \end{tabular}
# \end{table}

# For LGG
outDir <- paste(rootDir, "/Data", sep = "")
load(file = paste(outDir, "/LGG/LGG_lasso_min.Rdata", sep = "")) # return lasso_min
head(lasso_min)
nrow(lasso_min)
# > head(lasso_min)
# Active.Index Active.Coefficients lasso_gene
# 1            1         -0.12731255        ALB
# 2            3          0.25492445      CASP3
# 3            4          0.17735078       CDK1
# 4            5          0.15464554         CP
# 5            6         -0.04102212     CYP1A1
# 6            8         -0.03404832     MT-CO1
# > nrow(lasso_min)
# [1] 6

t <- lasso_min
t <- t[,c(3,2)]
colnames(t) <- c("Prognostic gene", "Coefficient")
head(t)
# > head(t)
# Prognostic gene Coefficient
# 1             ALB -0.12731255
# 2           CASP3  0.25492445
# 3            CDK1  0.17735078
# 4              CP  0.15464554
# 5          CYP1A1 -0.04102212
# 6          MT-CO1 -0.03404832

t[,2]<-format(as.numeric(t[,2]),big.mark=",")
print(latex.table.by(t, format.args = list(digits = 2, format = c("s","d"))), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)
# % latex table generated in R 4.0.2 by xtable 1.8-4 package
# % Tue Oct 11 15:33:10 2022
# \begin{table}[ht]
# \centering
# \begin{tabular}{|c|r|}
# \hline
# Prognostic gene & Coefficient \\
# \hline
# \cline{1-2}\multirow{ 1 }{*}{ ALB } & -0.12731255 \\
# \cline{1-2}\multirow{ 1 }{*}{ CASP3 } &  0.25492445 \\
# \cline{1-2}\multirow{ 1 }{*}{ CDK1 } &  0.17735078 \\
# \cline{1-2}\multirow{ 1 }{*}{ CP } &  0.15464554 \\
# \cline{1-2}\multirow{ 1 }{*}{ CYP1A1 } & -0.04102212 \\
# \cline{1-2}\multirow{ 1 }{*}{ MT-CO1 } & -0.03404832 \\
# \hline
# \end{tabular}
# \end{table}

#================================================================

#================================================================
# (13) Identify cancer types for tumor classification
# Katana
#================================================================

# 18 genes

# Get group information
outDir <- paste(rootDir, "/Data/Output", sep = "")
f18 <- paste(outDir, "/sur18.Rdata", sep = "")
load(file=f18) # return result1
group1=result1$group
distanceMatrix1=result1$distanceMatrix
head(group1)
length(group1)
table(group1)
distanceMatrix1[1:5,1:4]
nrow(distanceMatrix1)
ncol(distanceMatrix1)
# > head(group1)
# [1] 1 1 1 1 1 1
# > length(group1)
# [1] 6727
# > table(group1)
# group1
# 1    2
# 6537  190
# > distanceMatrix1[1:5,1:4]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01
# TCGA.OR.A5J1.01    2.775666e-01    1.667835e-04    4.458156e-05    2.478178e-04
# TCGA.OR.A5J2.01    1.667835e-04    2.775666e-01    1.285496e-05    5.112811e-04
# TCGA.OR.A5J3.01    4.458156e-05    1.285496e-05    2.775666e-01    1.841660e-05
# TCGA.OR.A5J5.01    2.478178e-04    5.112811e-04    1.841660e-05    2.775666e-01
# TCGA.OR.A5J6.01    5.714425e-04    1.546315e-05    1.094152e-04    6.690407e-05
# > nrow(distanceMatrix1)
# [1] 6727
# > ncol(distanceMatrix1)
# [1] 6727

# Get cancer types
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 59

# Merge the data
identical(rownames(distanceMatrix1), rownames(surdata))
# > identical(rownames(distanceMatrix1), rownames(surdata))
# [1] TRUE
patients <- surdata
patients$group <- group1
patients <- patients[, c("cancertype", "group")]
head(patients)
nrow(patients)
table(patients)
# > head(patients)
# cancertype group
# TCGA.OR.A5J1.01        ACC     1
# TCGA.OR.A5J2.01        ACC     1
# TCGA.OR.A5J3.01        ACC     1
# TCGA.OR.A5J5.01        ACC     1
# TCGA.OR.A5J6.01        ACC     1
# TCGA.OR.A5J7.01        ACC     1
# > nrow(patients)
# [1] 6727
# > table(patients)
# group
# cancertype    1    2
# ACC    76    0
# AML     0  172
# BLCA  144    0
# BRCA 1211    0
# COAD  274    3
# ESCA  179    0
# GBM   163    2
# KICH   91    0
# KIRC  507    1
# KIRP  188    1
# LGG   348    0
# LIHC  239    0
# LUAD  273    1
# LUSC  192    0
# OVCA  418    0
# PAAD  153    0
# PCPG  179    0
# PRAD  373    0
# SKCM  323    8
# STAD  408    1
# THCA  560    0
# UCEC  181    1
# UCS    57    0

# 12 genes

# Get group information
outDir <- paste(rootDir, "/Data/Output", sep = "")
f12 <- paste(outDir, "/sur12.Rdata", sep = "")
load(file=f12) # return result1
group1=result1$group
distanceMatrix1=result1$distanceMatrix
head(group1)
length(group1)
table(group1)
distanceMatrix1[1:5,1:4]
nrow(distanceMatrix1)
ncol(distanceMatrix1)
# > head(group1)
# [1] 1 1 1 1 1 1
# > length(group1)
# [1] 6727
# > table(group1)
# group1
# 1    2
# 6553  174
# > distanceMatrix1[1:5,1:4]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01
# TCGA.OR.A5J1.01    1.995860e-01    1.126702e-04    1.863774e-03    8.866152e-03
# TCGA.OR.A5J2.01    1.126702e-04    1.995860e-01    3.509226e-04    2.761649e-04
# TCGA.OR.A5J3.01    1.863774e-03    3.509226e-04    1.995860e-01    1.941792e-02
# TCGA.OR.A5J5.01    8.866152e-03    2.761649e-04    1.941792e-02    1.995860e-01
# TCGA.OR.A5J6.01    3.831460e-07    1.807445e-05    3.460988e-06    1.950237e-06
# > nrow(distanceMatrix1)
# [1] 6727
# > ncol(distanceMatrix1)
# [1] 6727

# Get cancer type
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 59

# Merge the data
identical(rownames(distanceMatrix1), rownames(surdata))
# > identical(rownames(distanceMatrix1), rownames(surdata))
# [1] TRUE
patients <- surdata
patients$group <- group1
patients <- patients[, c("cancertype", "group")]
head(patients)
nrow(patients)
table(patients)
# > head(patients)
# cancertype group
# TCGA.OR.A5J1.01        ACC     1
# TCGA.OR.A5J2.01        ACC     1
# TCGA.OR.A5J3.01        ACC     1
# TCGA.OR.A5J5.01        ACC     1
# TCGA.OR.A5J6.01        ACC     1
# TCGA.OR.A5J7.01        ACC     1
# > nrow(patients)
# [1] 6727
# > table(patients)
# group
# cancertype    1    2
# ACC    76    0
# AML   172    0
# BLCA  144    0
# BRCA 1211    0
# COAD  277    0
# ESCA  179    0
# GBM   165    0
# KICH   91    0
# KIRC  508    0
# KIRP  189    0
# LGG   348    0
# LIHC  239    0
# LUAD  274    0
# LUSC  192    0
# OVCA  418    0
# PAAD  152    1
# PCPG    7  172
# PRAD  373    0
# SKCM  331    0
# STAD  409    0
# THCA  560    0
# UCEC  181    1
# UCS    57    0

#================================================================

#================================================================
# (14) Latex table for genes up down expressed in tumour
# Run in local
#================================================================

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
# > y[1:5,1:4]
# gene frequency ACC  AML
# 1     CDK1        18  UP   UP
# 2      ALB        13   0 DOWN
# 3    AP1S1        11   0 DOWN
# 4    CASP3         9   0    0
# 5 MAP1LC3A         8   0 DOWN
# > nrow(y)
# [1] 56
# > ncol(y)
# [1] 27

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# > head(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# > nrow(uni_gene)
# [1] 35

# Get necessary data
dat <- y
dat <- dat[which(dat$gene %in% uni_gene$x),]
dat <- dat[,c(1,2,26,27)]
dat$perUp <- paste(round((dat$numUp/dat$frequency)*100,1), "\\%", sep = "")
dat$perDown <- paste(round((dat$numDown/dat$frequency)*100,1), "\\%", sep = "")
dat$group <- ifelse(dat$numUp > dat$numDown, "Up", ifelse(dat$numUp < dat$numDown, "Down", "Varied"))
dat <- dat[,c(1,2,5,6,7)]
colnames(dat) <- c('CCG', 'Number of cancer types', 'Percentage for up expressed', 'Percentage for down expressed', 'Change')
head(dat)
# > head(dat)
# CCG Number of cancer types Percentage for up expressed
# 1     CDK1                     18                        100%
# 3    AP1S1                     11                       90.9%
# 4    CASP3                      9                        100%
# 5 MAP1LC3A                      8                       12.5%
# 6     SNCA                      8                         25%
# 7  TMPRSS6                      7                       57.1%
# Percentage for down expressed Change
# 1                            0%     Up
# 3                          9.1%     Up
# 4                            0%     Up
# 5                         87.5%   Down
# 6                           75%   Down
# 7                         42.9%     Up

# Latex table
dat[,2]<-format(as.numeric(dat[,2]),big.mark=",")
print(latex.table.by(dat, format.args = list(digits = 2, format = c("s", "d", "s", "s", "s"))), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

#================================================================

#================================================================
# (15) SNV
# Local
#================================================================

# Read critical copper genes
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# All 56 critical copper genes
geneList56 <- x$gene

# 39 popular critical cooper genes
geneList39 <- x[which(x$frequency > 1),]$gene

# 20 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList20 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# > head(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# > nrow(uni_gene)
# [1] 35

# Cox genes
geneList18 <- geneList20[which(geneList20 %in% uni_gene$x)]
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]
print("35 genes")
print(paste(uni_gene$x, collapse = ", "))
print("18 genes")
print(paste(geneList18, collapse = ", "))
print("12 genes")
print(paste(geneList12, collapse = ", "))
# > print("35 genes")
# [1] "35 genes"
# > print(paste(uni_gene$x, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, MAP1LC3A, SNCA, TMPRSS6, MAPT, GSK3B, JUN, APP, CYP1A1, COX17, XIAP, FOXO1, ARF1, GPC1, AOC3, SORD, PRNP, ATP7A, SP1, MT-CO1, DBH, SLC11A2, ANG, S100A12, ATP6AP1, ADAM10, MT-CO2, CP, PRND, AQP1, MT1X, IL1A, XAF1"
# > print("18 genes")
# [1] "18 genes"
# > print(paste(geneList18, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, ARF1, GPC1, SORD, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"
# > print("12 genes")
# [1] "12 genes"
# > print(paste(geneList12, collapse = ", "))
# [1] "MAP1LC3A, SNCA, MAPT, JUN, CYP1A1, AOC3, PRNP, DBH, S100A12, AQP1, MT1X, XAF1"

#================================================================

#================================================================
# (16) Copper genes - Genes related to copper metabolism
# Local
#================================================================
# https://bioconductor.org/packages/release/data/experiment/vignettes/msigdb/inst/doc/msigdb.html

library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(fgsea)

# Local
rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# rootDir="/srv/scratch/z3538133/002NetworkAnalysis" # And put the input files in "rootDir/Data"
setwd(rootDir)
source(paste(rootDir, "/Script/ProposedMethod_Functions.R", sep=""))

# Pathways
pathways <- c("WP_COPPER_HOMEOSTASIS", "HP_DECREASED_CIRCULATING_COPPER_CONCENTRATION",
              "HP_ABNORMAL_CIRCULATING_COPPER_CONCENTRATION", "GOMF_COPPER_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
              "GOMF_COPPER_ION_BINDING", "GOMF_COPPER_CHAPERONE_ACTIVITY",
              "GOBP_RESPONSE_TO_COPPER_ION", "GOBP_DETOXIFICATION_OF_COPPER_ION",
              "GOBP_COPPER_ION_TRANSPORT", "GOBP_COPPER_ION_TRANSMEMBRANE_TRANSPORT",
              "GOBP_COPPER_ION_IMPORT", "GOBP_COPPER_ION_HOMEOSTASIS",
              "GOBP_CELLULAR_RESPONSE_TO_COPPER_ION", "GOBP_CELLULAR_COPPER_ION_HOMEOSTASIS")

GMTdirLoc <- paste(rootDir, "/Data/msigdb_v2022.1.Hs_files_to_download_locally/msigdb_v2022.1.Hs_GMTs/", sep = "") # Directory containing all MSigDB .gmt files
files <- c("c2.cp.wikipathways.v2022.1.Hs.symbols.gmt", "c5.hpo.v2022.1.Hs.symbols.gmt", "c5.go.mf.v2022.1.Hs.symbols.gmt", "c5.go.bp.v2022.1.Hs.symbols.gmt")

r <- matrix(NA, nrow = 14, ncol = 2)
ind <- 0
allGenes <- c()
for (i in 1:4) {
  pathwaysSub <- gmtPathways(paste(GMTdirLoc, files[i], sep = ""))
  ids <- which(names(pathwaysSub) %in% pathways)
  for (j in 1:length(ids)) {
    ind <- ind + 1
    r[ind,1] <- names(pathwaysSub)[ids[j]]
    r[ind,2] <- paste(pathwaysSub[[ids[j]]], collapse=", ")
    allGenes <- c(allGenes, pathwaysSub[[ids[j]]])
  }
}

# Test copper genes
# Test pathways
identical(pathways[order(pathways, decreasing = TRUE)], r[order(r[,1], decreasing = TRUE),1])
# [1] TRUE
# test 133 genes
# copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'CA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')
copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'COA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')
allGenes <- unique(allGenes)
# # Gene CA6 in copper.genelist should be COA6
# copper.genelist <- copper.genelist[-which(copper.genelist[] == "CA6")]
# copper.genelist <- c(copper.genelist, "COA6")
copper.genelist <- copper.genelist[order(copper.genelist, decreasing = TRUE)]
allGenes <- allGenes[order(allGenes, decreasing = TRUE)]
identical(copper.genelist,allGenes)
# [1] TRUE

colnames(r) <- c("Pathway", "Genes")
write.xlsx(as.data.frame(r), paste(rootDir, "/Data/Output/CopperMetabolismPathways.xlsx", sep = ""))
write.xlsx(as.data.frame(copper.genelist), paste(rootDir, "/Data/Output/CopperMetabolismGenes.xlsx", sep = ""))

#================================================================
# Pathways from Toni
# Local
rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# rootDir="/srv/scratch/z3538133/002NetworkAnalysis" # And put the input files in "rootDir/Data"
setwd(rootDir)
source(paste(rootDir, "/Script/ProposedMethod_Functions.R", sep=""))

# Pathways
pathways <- c("WP_COPPER_HOMEOSTASIS", "HP_DECREASED_CIRCULATING_COPPER_CONCENTRATION",
              "HP_ABNORMAL_CIRCULATING_COPPER_CONCENTRATION", "GOMF_COPPER_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY",
              "GOMF_COPPER_ION_BINDING", "GOMF_COPPER_CHAPERONE_ACTIVITY",
              "GOBP_RESPONSE_TO_COPPER_ION", "GOBP_DETOXIFICATION_OF_COPPER_ION",
              "GOBP_COPPER_ION_TRANSPORT", "GOBP_COPPER_ION_TRANSMEMBRANE_TRANSPORT",
              "GOBP_COPPER_ION_IMPORT", "GOBP_COPPER_ION_HOMEOSTASIS",
              "GOBP_CELLULAR_RESPONSE_TO_COPPER_ION", "GOBP_CELLULAR_COPPER_ION_HOMEOSTASIS")

GMTdirLoc <- paste(rootDir, "/Data/copper_genesets 2/", sep = "") # Directory containing all MSigDB .gmt files

r <- matrix(NA, nrow = 14, ncol = 2)
ind <- 0
allGenes <- c()
for (i in 1:14) {
  pathwaysSub <- gmtPathways(paste(GMTdirLoc, pathways[i], ".v7.5.1.gmt", sep = ""))
  ind <- ind + 1
  r[ind,1] <- names(pathwaysSub)[1]
  r[ind,2] <- paste(pathwaysSub[[1]], collapse=", ")
  allGenes <- c(allGenes, pathwaysSub[[1]])
}

# Test copper genes
# Test pathways
identical(pathways[order(pathways, decreasing = TRUE)], r[order(r[,1], decreasing = TRUE),1])
# [1] TRUE
# test 133 genes
# copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'CA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')
copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'COA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')
allGenes <- unique(allGenes)
# # Gene CA6 in copper.genelist should be COA6
# copper.genelist <- copper.genelist[-which(copper.genelist[] == "CA6")]
# copper.genelist <- c(copper.genelist, "COA6")
copper.genelist <- copper.genelist[order(copper.genelist, decreasing = TRUE)]
allGenes <- allGenes[order(allGenes, decreasing = TRUE)]
identical(copper.genelist,allGenes)
# [1] TRUE

colnames(r) <- c("Pathway", "Genes")
write.xlsx(as.data.frame(r), paste(rootDir, "/Data/Output/CopperMetabolismPathways.xlsx", sep = ""))
write.xlsx(as.data.frame(copper.genelist), paste(rootDir, "/Data/Output/CopperMetabolismGenes.xlsx", sep = ""))

#================================================================

#================================================================
# (17) Expression of gene SLC31A1
# Katana
#================================================================
# Katana
# rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"
rootDir="/srv/scratch/z3538133/002NetworkAnalysis" # And put the input files in "rootDir/Data"
setwd(rootDir)
source(paste(rootDir, "/Script/ProposedMethod_Functions.R", sep=""))

setwd(paste(rootDir, "/Data/SNV", sep = ""))

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

SNV_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del",
               "In_Frame_Ins")

set35Genes <- c("ATP7A", "CP", "TMPRSS6", "MAPT", "DBH",
                "APP", "AOC3", "ADAM10", "SP1", "CYP1A1",
                "XIAP", "FOXO1", "ATP6AP1", "SLC11A2", "GSK3B",
                "GPC1", "PRNP", "AP1S1", "AQP1", "JUN",
                "CDK1", "ARF1", "CASP3", "SORD", "IL1A",
                "MAP1LC3A", "XAF1", "PRND", "SNCA", "S100A12",
                "ANG", "COX17", "MT1X") # 33 left, 2 genes, MT-CO1 and MT-CO2, were removed due to the missing data

top10Genes <- c("ATP7A", "CP", "APP", "MAPT", "DBH",
                "TMPRSS6", "AOC3", "ADAM10", "SP1", "XIAP")

geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= top10Genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
# > geneIDs
# SYMBOL          GENEID
# 1    ATP7A ENSG00000165240
# 2       CP ENSG00000047457
# 3      APP ENSG00000142192
# 4     MAPT ENSG00000186868
# 5     MAPT ENSG00000276155
# 6     MAPT ENSG00000277956
# 7      DBH ENSG00000123454
# 8  TMPRSS6 ENSG00000187045
# 9     AOC3 ENSG00000131471
# 10  ADAM10 ENSG00000137845
# 11     SP1 ENSG00000185591
# 12    XIAP ENSG00000101966
# 13    XIAP          LRG_19

geneID35s <- ensembldb::select(EnsDb.Hsapiens.v79, keys= set35Genes, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

n <- length(subtypes)
for (i in 1:n) {
  cancertype <- subtypes[i]
  # Read data
  dat <- as.data.frame(fread(paste(cancertype, "_maf_data.IdTrans.tsv", sep = "")))
  # Only get necessary data
  dat <- dat[,c("Variant_Classification", "Tumor_Sample_Barcode", "Gene")]
  # > head(dat)
  # Variant_Classification         Tumor_Sample_Barcode            Gene
  # 1      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000122375
  # 2                 Intron TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000090615
  # 3      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000111796
  # 4      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000165821
  # 5                 Silent TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000137841
  # 6        Frame_Shift_Del TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000169758
  dat <- dat[which(dat$Variant_Classification %in% SNV_types),]
  dat$CancerType <- cancertype
  if(i == 1) {
    r <- dat
  } else {
    r <- rbind(r,dat)
  }
}

# Remove items from geneIDs which are not in r
idx <- which(!(geneIDs$GENEID %in% r$Gene))
top10GeneIDs <- geneIDs[-idx,]
# > top10GeneIDs
# SYMBOL          GENEID
# 1    ATP7A ENSG00000165240
# 2       CP ENSG00000047457
# 3      APP ENSG00000142192
# 4     MAPT ENSG00000186868
# 7      DBH ENSG00000123454
# 8  TMPRSS6 ENSG00000187045
# 9     AOC3 ENSG00000131471
# 10  ADAM10 ENSG00000137845
# 11     SP1 ENSG00000185591
# 12    XIAP ENSG00000101966
# top10Genes <- c("ATP7A", "CP", "APP", "MAPT", "DBH",
#                 "TMPRSS6", "AOC3", "ADAM10", "SP1", "XIAP")

idx <- which(!(geneID35s$GENEID %in% r$Gene))
top35GeneIDs <- geneID35s[-idx,]
head(top35GeneIDs)
nrow(top35GeneIDs)
# > head(top35GeneIDs)
# SYMBOL          GENEID
# 1   ATP7A ENSG00000165240
# 2      CP ENSG00000047457
# 3 TMPRSS6 ENSG00000187045
# 4    MAPT ENSG00000186868
# 7     DBH ENSG00000123454
# 8     APP ENSG00000142192
# > nrow(top35GeneIDs)
# [1] 33

# Mutations for all genes
saveRDS(r, "SNVAllGenes.rds")
# r <- readRDS("SNVAllGenes.rds")

# 35 genes
saveRDS(r[which(r$Gene %in% top35GeneIDs$GENEID),], "SNV35Genes.rds")

# Mutations for 10 genes
saveRDS(r[which(r$Gene %in% top10GeneIDs$GENEID),], "SNV10Genes.rds")

# Mutations for ATP7A
saveRDS(r[which(r$Gene == "ENSG00000165240"),], "SNVATP7A.rds")

SNVAllGenes <- readRDS("SNVAllGenes.rds")
head(SNVAllGenes)
# > head(SNVAllGenes)
# Variant_Classification         Tumor_Sample_Barcode            Gene
# 1      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000122375
# 3      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000111796
# 4      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000165821
# 6        Frame_Shift_Del TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000169758
# 7      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000104731
# 9      Missense_Mutation TCGA-OR-A5J1-01A-11D-A29I-10 ENSG00000130935
# CancerType
# 1        ACC
# 3        ACC
# 4        ACC
# 6        ACC
# 7        ACC
# 9        ACC
SNV35Genes <- readRDS("SNV35Genes.rds")
head(SNV35Genes)
SNV10Genes <- readRDS("SNV10Genes.rds")
head(SNV10Genes)
# > head(SNV10Genes)
# Variant_Classification         Tumor_Sample_Barcode            Gene
# 452        Missense_Mutation TCGA-OR-A5J5-01A-11D-A29I-10 ENSG00000131471
# 606          Frame_Shift_Del TCGA-OR-A5J5-01A-11D-A29I-10 ENSG00000187045
# 1030       Missense_Mutation TCGA-OR-A5J8-01A-11D-A29I-10 ENSG00000123454
# 9226       Missense_Mutation TCGA-OR-A5LL-01A-11D-A29I-10 ENSG00000047457
# 20912      Missense_Mutation TCGA-4Z-AA7M-01A-11D-A391-08 ENSG00000137845
# 6059            In_Frame_Del TCGA-4Z-AA84-01A-11D-A391-08 ENSG00000185591
# CancerType
# 452          ACC
# 606          ACC
# 1030         ACC
# 9226         ACC
# 20912       BLCA
# 6059        BLCA
SNVATP7A <- readRDS("SNVATP7A.rds")
head(SNVATP7A)
# > head(SNVATP7A)
# Variant_Classification         Tumor_Sample_Barcode            Gene
# 109072      Missense_Mutation TCGA-BT-A20J-01A-11D-A14W-08 ENSG00000165240
# 62121       Missense_Mutation TCGA-FD-A3N6-01A-11D-A21A-08 ENSG00000165240
# 64539       Missense_Mutation TCGA-FD-A3SQ-01A-21D-A22Z-08 ENSG00000165240
# 89908       Missense_Mutation TCGA-GU-A763-01A-11D-A32B-08 ENSG00000165240
# 118677      Missense_Mutation TCGA-UY-A78L-01A-12D-A339-08 ENSG00000165240
# 133255      Missense_Mutation TCGA-XF-A9T6-01A-11D-A42E-08 ENSG00000165240
# CancerType
# 109072       BLCA
# 62121        BLCA
# 64539        BLCA
# 89908        BLCA
# 118677       BLCA
# 133255       BLCA

# For all genes, ~ 20K genes
# SNVAllGenes
# Number of patients
npatientsAll <- length(unique(SNVAllGenes$Tumor_Sample_Barcode))
npatientsAll
# > npatientsAll
# [1] 8523

# 35 genes
npatients35 <- length(unique(SNV35Genes$Tumor_Sample_Barcode))
npatients35
# > npatients35
# [1] 1280

# For 10 genes
# SNV10Genes
# Number of patients
npatients10 <- length(unique(SNV10Genes$Tumor_Sample_Barcode))
npatients10
# > npatients10
# [1] 887

# For SNVATP7A
# SNVATP7A
# Number of patients
npatientsATP7A <- length(unique(SNVATP7A$Tumor_Sample_Barcode))
npatientsATP7A
# > npatientsATP7A
# [1] 202

# Load the gene expression data
loadData=function(gene) {
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- unique(subtypes$cancer.type)
  subtypes <- c(subtypes, "PAAD")
  subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
  subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
  n <- length(subtypes)
  Cancer_type <- c()
  Type <- c() # Tumour or Normal
  Expression <- c()
  Sample <- c() # Sample id
  for (i in 1:n) {
    cancertype <- subtypes[i]
    outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
    load(paste(outDir, "/gtex.", cancertype, ".raw.counts.RData", sep = "")) # return gtex.PAN.raw.counts
    # transform
    dat <- gtex.PAN.raw.counts
    rownames(dat) <- dat$gene
    dat <- dat[,-1]
    dat <- dat + 1
    dat <- log(dat, base = 2)
    iGene <- which(rownames(dat) == gene)
    r <- str_detect(colnames(dat), "TCGA")
    nTumor <- length(which(r == TRUE))
    nNormal <- ncol(dat) - nTumor
    
    Cancer_type <- c(Cancer_type, rep(cancertype, nTumor + nNormal))
    Type <- c(Type, rep(c("Tumour", "Normal"), times = c(nTumor, nNormal)))
    Expression <- c(Expression, as.numeric(dat[iGene,]))
    Sample <- c(Sample, colnames(dat))
  }
  data=data.frame(Cancer_type, Type ,  Expression, Sample)
  
  return(data)
}

data <- loadData("SLC31A1")
head(data)
nrow(data)
ncol(data)
# > head(data)
# Cancer_type   Type Expression          Sample
# 1         ACC Tumour     9.1898 TCGA.OR.A5K3.01
# 2         ACC Tumour    11.5430 TCGA.OR.A5J2.01
# 3         ACC Tumour     8.8887 TCGA.OR.A5LN.01
# 4         ACC Tumour    12.6425 TCGA.OR.A5KY.01
# 5         ACC Tumour    10.5038 TCGA.OR.A5LG.01
# 6         ACC Tumour    10.0334 TCGA.OR.A5JC.01
# > nrow(data)
# [1] 11329
# > ncol(data)
# [1] 4

# Get tumor only
data2 <- data[which(data$Type == "Tumour"),]
head(data2)
nrow(data2)
ncol(data2)
# > head(data2)
# Cancer_type   Type Expression          Sample
# 1         ACC Tumour     9.1898 TCGA.OR.A5K3.01
# 2         ACC Tumour    11.5430 TCGA.OR.A5J2.01
# 3         ACC Tumour     8.8887 TCGA.OR.A5LN.01
# 4         ACC Tumour    12.6425 TCGA.OR.A5KY.01
# 5         ACC Tumour    10.5038 TCGA.OR.A5LG.01
# 6         ACC Tumour    10.0334 TCGA.OR.A5JC.01
# > nrow(data2)
# [1] 6732
# > ncol(data2)
# [1] 4

# Remove samples with ATP7A mutation from SNV10Genes and SNVAllGenes
patientsATP7A <- unique(SNVATP7A$Tumor_Sample_Barcode)
SNV10Genes <- SNV10Genes[which(!(SNV10Genes$Tumor_Sample_Barcode %in% patientsATP7A)),]
length(unique(SNV10Genes$Tumor_Sample_Barcode))
# > length(unique(SNV10Genes$Tumor_Sample_Barcode))
# [1] 685
SNVAllGenes <- SNVAllGenes[which(!(SNVAllGenes$Tumor_Sample_Barcode %in% patientsATP7A)),]
length(unique(SNVAllGenes$Tumor_Sample_Barcode))
# > length(unique(SNVAllGenes$Tumor_Sample_Barcode))
# [1] 8321

# Get only patients with ATP7A mutation
# For SNVATP7A
patientsATP7A <- unique(SNVATP7A$Tumor_Sample_Barcode)
head(patientsATP7A)
length(patientsATP7A)
# > head(patientsATP7A)
# [1] "TCGA-BT-A20J-01A-11D-A14W-08" "TCGA-FD-A3N6-01A-11D-A21A-08"
# [3] "TCGA-FD-A3SQ-01A-21D-A22Z-08" "TCGA-GU-A763-01A-11D-A32B-08"
# [5] "TCGA-UY-A78L-01A-12D-A339-08" "TCGA-XF-A9T6-01A-11D-A42E-08"
# > length(patientsATP7A)
# [1] 202
# Format sample id
patientsATP7A <- substr(patientsATP7A, 1, 15)
patientsATP7A <- gsub('-','.', patientsATP7A)
data2ATP7A <- data2[which(data2$Sample %in% patientsATP7A),]
head(data2ATP7A)
nrow(data2ATP7A)
summary(data2ATP7A$Expression)
# > head(data2ATP7A)
# Cancer_type   Type Expression          Sample
# 860         BLCA Tumour    10.7322 TCGA.FD.A3SQ.01
# 949         BLCA Tumour    11.2064 TCGA.BT.A20J.01
# 1019        BRCA Tumour    10.7108 TCGA.A2.A0CR.01
# 1036        BRCA Tumour    13.1542 TCGA.B6.A0WZ.01
# 1067        BRCA Tumour    11.8321 TCGA.S3.AA14.01
# 1122        BRCA Tumour    13.8388 TCGA.D8.A27V.01
# > nrow(data2ATP7A)
# [1] 119
# > summary(data2ATP7A$Expression)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 9.867  11.164  11.733  11.770  12.386  13.938

# Get patients without ATP7A mutation among 10 genes
patients10Genes <- unique(SNV10Genes$Tumor_Sample_Barcode)
head(patients10Genes)
length(patients10Genes)
# > head(patients10Genes)
# [1] "TCGA-OR-A5J5-01A-11D-A29I-10" "TCGA-OR-A5J8-01A-11D-A29I-10"
# [3] "TCGA-OR-A5LL-01A-11D-A29I-10" "TCGA-4Z-AA7M-01A-11D-A391-08"
# [5] "TCGA-4Z-AA84-01A-11D-A391-08" "TCGA-BL-A5ZZ-01A-31D-A30E-08"
# > length(patients10Genes)
# [1] 685
# Format sample id
patients10Genes <- substr(patients10Genes, 1, 15)
patients10Genes <- gsub('-','.', patients10Genes)
data210genes <- data2[which(data2$Sample %in% patients10Genes),]
head(data210genes)
nrow(data210genes)
summary(data210genes$Expression)
# > head(data210genes)
# Cancer_type   Type Expression          Sample
# 16          ACC Tumour    11.1485 TCGA.OR.A5LL.01
# 38          ACC Tumour    10.2180 TCGA.OR.A5J5.01
# 48          ACC Tumour    12.2009 TCGA.OR.A5J8.01
# 836        BLCA Tumour    11.4757 TCGA.G2.A2EO.01
# 843        BLCA Tumour    12.2949 TCGA.GV.A3QI.01
# 865        BLCA Tumour    12.6432 TCGA.GV.A3JX.01
# > nrow(data210genes)
# [1] 402
# > summary(data210genes$Expression)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 8.234  11.225  11.807  11.802  12.343  14.925

# Get patients without ATP7A mutation in all genes
patientsAll <- unique(SNVAllGenes$Tumor_Sample_Barcode)
head(patientsAll)
length(patientsAll)
# > head(patientsAll)
# [1] "TCGA-OR-A5J1-01A-11D-A29I-10" "TCGA-OR-A5J2-01A-11D-A29I-10"
# [3] "TCGA-OR-A5J3-01A-11D-A29I-10" "TCGA-OR-A5J4-01A-11D-A29I-10"
# [5] "TCGA-OR-A5J5-01A-11D-A29I-10" "TCGA-OR-A5J6-01A-31D-A29I-10"
# > length(patientsAll)
# [1] 8321
# Format sample id
patientsAll <- substr(patientsAll, 1, 15)
patientsAll <- gsub('-','.', patientsAll)
data2All <- data2[which(data2$Sample %in% patientsAll),]
head(data2All)
nrow(data2All)
summary(data2All$Expression)
# > head(data2All)
# Cancer_type   Type Expression          Sample
# 1         ACC Tumour     9.1898 TCGA.OR.A5K3.01
# 2         ACC Tumour    11.5430 TCGA.OR.A5J2.01
# 3         ACC Tumour     8.8887 TCGA.OR.A5LN.01
# 4         ACC Tumour    12.6425 TCGA.OR.A5KY.01
# 5         ACC Tumour    10.5038 TCGA.OR.A5LG.01
# 6         ACC Tumour    10.0334 TCGA.OR.A5JC.01
# > nrow(data2All)
# [1] 5401
# > summary(data2All$Expression)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 5.358  10.903  11.576  11.545  12.210  15.121

f <- paste(rootDir, "/Data/Output/SLC31A1.pdf", sep = "")
pdf(file = f, width = 5, height =  5)

data2ATP7A$Mutated <- "MutatedATP7A"
data2All$Mutated <- "NotMutatedATP7A"
combinedData <- rbind(data2ATP7A, data2All)
head(combinedData)
nrow(combinedData)
ncol(combinedData)
# > head(combinedData)
# Cancer_type   Type Expression          Sample      Mutated
# 860         BLCA Tumour    10.7322 TCGA.FD.A3SQ.01 MutatedATP7A
# 949         BLCA Tumour    11.2064 TCGA.BT.A20J.01 MutatedATP7A
# 1019        BRCA Tumour    10.7108 TCGA.A2.A0CR.01 MutatedATP7A
# 1036        BRCA Tumour    13.1542 TCGA.B6.A0WZ.01 MutatedATP7A
# 1067        BRCA Tumour    11.8321 TCGA.S3.AA14.01 MutatedATP7A
# 1122        BRCA Tumour    13.8388 TCGA.D8.A27V.01 MutatedATP7A
# > nrow(combinedData)
# [1] 5520
# > ncol(combinedData)
# [1] 5

boxplot(Expression~Mutated,data=combinedData, main="",
        xlab="Mutated in ATP7A", ylab="SLC31A1 Expression")

dev.off()

#================================================================

#================================================================
# (18) Survival analysis, using 18 genes and 12 genes for each cancer type
# Run in server
#================================================================

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# We do not use lasso for pan-cancer, only use univariate Cox
# # Read lasso genes
# load(file = paste(outDir, "/lasso_min.Rdata", sep = "")) # return lasso_min
# nrow(lasso_min)
# head(lasso_min)
# # > nrow(lasso_min)
# # [1] 32
# # > head(lasso_min)
# # Active.Index Active.Coefficients lasso_gene
# # 1            1          0.05136781     ADAM10
# # 2            2          0.12289435        ANG
# # 3            3         -0.02669716       AOC3
# # 4            4         -0.01338193      AP1S1
# # 5            5         -0.16595332        APP
# # 6            6          0.02119514       AQP1

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# > head(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# > nrow(uni_gene)
# [1] 35

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
# > geneList18
# [1] "CDK1"    "AP1S1"   "CASP3"   "TMPRSS6" "GSK3B"   "APP"     "COX17"
# [8] "XIAP"    "ARF1"    "GPC1"    "SORD"    "ATP7A"   "SP1"     "MT-CO1"
# [15] "SLC11A2" "ATP6AP1" "ADAM10"  "CP"
# f18 <- paste(outDir, "/sur18.pdf", sep = "")
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]
# > geneList12
# [1] "MAP1LC3A" "SNCA"     "MAPT"     "JUN"      "CYP1A1"   "AOC3"
# [7] "PRNP"     "DBH"      "S100A12"  "AQP1"     "MT1X"     "XAF1"
# f12 <- paste(outDir, "/sur12.pdf", sep = "")
# K	Number of nearest neighbors, default is 20
# surAnalysis(surdata, geneList18, numGroup=2, K=20, alpha=0.5, f18, w = 9, h = 6)
# surAnalysis(surdata, geneList12, numGroup=2, K=20, alpha=0.5, f12, w = 9, h = 6)
print(paste(geneList18, collapse = ", "))
print(paste(geneList12, collapse = ", "))
# > print(paste(geneList18, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, ARF1, GPC1, SORD, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"
# > print(paste(geneList12, collapse = ", "))
# [1] "MAP1LC3A, SNCA, MAPT, JUN, CYP1A1, AOC3, PRNP, DBH, S100A12, AQP1, MT1X, XAF1"

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
n <- length(subtypes)

subDir <- "surCancerType1812"
mainDir <- paste(rootDir, "/Data/Output", sep = "")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

pList <- matrix(NA, nrow = n, ncol = 5)
colnames(pList) <- c("Cancer_type", "pvalue_18_genes", "pvalue_18_genes_significant", "pvalue_12_genes", "pvalue_12_genes_significant")

# pdf file
for (iType in 1:n) {
  cancertype <- subtypes[iType]
  outDir <- paste(rootDir, "/Data/Output/surCancerType1812", sep = "")
  dat <- surdata[which(surdata$cancertype == cancertype),]
  
  # # Critical copper genes
  # geneList <- criticalCopperGenes[which(criticalCopperGenes[,iType+2] == 1),]$gene
  f12 <- paste(outDir, "/", cancertype, "_sur12.pdf", sep = "")
  f18 <- paste(outDir, "/", cancertype, "_sur18.pdf", sep = "")
  
  p12 <- surAnalysis(dat, geneList12, numGroup=2, K=20, alpha=0.5, f12, w = 9, h = 6)
  p18 <- surAnalysis(dat, geneList18, numGroup=2, K=20, alpha=0.5, f18, w = 9, h = 6)
  
  pList[iType,1] <- cancertype
  pList[iType,2] <- p18
  pList[iType,3] <- if(p18 <= 0.05) {"Yes"} else {""}
  pList[iType,4] <- p12
  pList[iType,5] <- if(p12 <= 0.05) {"Yes"} else {""}
  
}

write.csv(pList, paste(rootDir, "/Data/Output/surCancerType1812/pList.csv", sep = ""))

#================================================================

#================================================================
# (19) Survival analysis, using critical copper genes for each cancer type
# Run in server
#================================================================

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
# > y[1:4,1:5]
# gene frequency ACC  AML BLCA
# 1  CDK1        18  UP   UP   UP
# 2   ALB        13   0 DOWN    0
# 3 AP1S1        11   0 DOWN    0
# 4 CASP3         9   0    0    0

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
n <- length(subtypes)

subDir <- "surCancerTypeUseCopperGenes"
mainDir <- paste(rootDir, "/Data/Output", sep = "")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

pList <- matrix(NA, nrow = n, ncol = 7)
colnames(pList) <- c("Cancer_type", "pvalue_up_genes", "pvalue_up_genes_significant", "pvalue_down_genes", "pvalue_down_genes_significant", "pvalue_all_genes", "pvalue_all_genes_significant")

# pdf file
for (iType in 1:n) {
  cancertype <- subtypes[iType]
  outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
  dat <- surdata[which(surdata$cancertype == cancertype),]
  
  # Critical copper genes
  criticalCopperGenes <- y
  upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,iType+2] == 'UP'),]$gene
  downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,iType+2] == 'DOWN'),]$gene
  geneList <- criticalCopperGenes[which(criticalCopperGenes[,iType+2] %in% c("UP", "DOWN")),]$gene
  fup <- paste(outDir, "/", cancertype, "_surUp.pdf", sep = "")
  fdown <- paste(outDir, "/", cancertype, "_surDown.pdf", sep = "")
  f <- paste(outDir, "/", cancertype, "_sur.pdf", sep = "")
  
  if(length(upGeneList) > 0) {
    pUp <- surAnalysis(dat, upGeneList, numGroup=2, K=20, alpha=0.5, fup, w = 9, h = 6)  
  } else {
    pUp <- NA
  }
  if(length(downGeneList) > 0) {
    pDown <- surAnalysis(dat, downGeneList, numGroup=2, K=20, alpha=0.5, fdown, w = 9, h = 6)
  } else {
    pDown <- NA
  }
  p <- surAnalysis(dat, geneList, numGroup=2, K=20, alpha=0.5, f, w = 9, h = 6)
  
  pList[iType,1] <- cancertype
  pList[iType,2] <- pUp
  if(!is.na(pUp)) {
    pList[iType,3] <- if(pUp <= 0.05) {"Yes"} else {""}  
  }
  pList[iType,4] <- pDown
  if(!is.na(pDown)) {
    pList[iType,5] <- if(pDown <= 0.05) {"Yes"} else {""}  
  }
  pList[iType,6] <- p
  pList[iType,7] <- if(p <= 0.05) {"Yes"} else {""}
  
}

write.csv(pList, paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes/pListUseCopperGenes.csv", sep = ""))

#----------------------------------------
# For triple negative breast cancer

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y[1:4,1:5]
# > y[1:4,1:5]
# gene frequency ACC  AML BLCA
# 1  CDK1        18  UP   UP   UP
# 2   ALB        13   0 DOWN    0
# 3 AP1S1        11   0 DOWN    0
# 4 CASP3         9   0    0    0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get only triple negative breast cancer (i.e., Basal)

# Get cancer sub types
subtypes <- c("BRCA")
all <- PanCancerAtlas_subtypes()
PAN.path.subtypes <- all[which(all$cancer.type %in% subtypes),]
PAN.path.subtypes[1:5,1:6]
nrow(PAN.path.subtypes)
# # A tibble: 5 × 6
# pan.samplesID                cancer.type Subtype_mRNA Subtyp…¹ Subty…² Subty…³
# <chr>                        <chr>       <chr>        <chr>      <int> <chr>
#   1 TCGA-E2-A158-11A-22R-A12D-07 BRCA        Normal       NA            NA NA
# 2 TCGA-BH-A0DD-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 3 TCGA-BH-A1EO-11A-31R-A137-07 BRCA        LumA         NA            NA NA
# 4 TCGA-BH-A0B5-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 5 TCGA-A7-A13G-11A-51R-A13Q-07 BRCA        LumA         NA            NA NA
# # … with abbreviated variable names ¹​Subtype_DNAmeth, ²​Subtype_protein,
# #   ³​Subtype_miRNA
# [1] 1218
PAN.subtypes <- PAN.path.subtypes
i <- sapply(PAN.subtypes, is.factor)
PAN.subtypes[i] <- lapply(PAN.subtypes[i], as.character)
table(PAN.subtypes$Subtype_mRNA)
# table(PAN.subtypes$Subtype_mRNA)
# 
# Basal   Her2   LumA   LumB Normal
# 193     82    581    219    143

# # make list of column IDs
# tumour.pt.ID <- as.list(PAN.subtypes$pan.samplesID)
# 
# # replace - with . to match the colnames in the counts dataframe
# tmp1 <- str_replace_all(tumour.pt.ID, "-", ".")

PAN.subtypes <- as.data.frame(PAN.subtypes)
PAN.subtypes$pan.samplesID <- substr(PAN.subtypes$pan.samplesID, 1, 15)
PAN.subtypes[1:5,1:6]
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA-E2-A158-11        BRCA       Normal            <NA>              NA
# 2 TCGA-BH-A0DD-11        BRCA         LumA            <NA>              NA
# 3 TCGA-BH-A1EO-11        BRCA         LumA            <NA>              NA
# 4 TCGA-BH-A0B5-11        BRCA         LumA            <NA>              NA
# 5 TCGA-A7-A13G-11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
PAN.subtypes$pan.samplesID <- str_replace_all(PAN.subtypes$pan.samplesID, "-", ".")
PAN.subtypes[1:5,1:6]
nrow(PAN.subtypes)
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA.E2.A158.11        BRCA       Normal            <NA>              NA
# 2 TCGA.BH.A0DD.11        BRCA         LumA            <NA>              NA
# 3 TCGA.BH.A1EO.11        BRCA         LumA            <NA>              NA
# 4 TCGA.BH.A0B5.11        BRCA         LumA            <NA>              NA
# 5 TCGA.A7.A13G.11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
#   [1] 1218

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subtype <- "BRCA"

subDir <- "surCancerTypeUseCopperGenes"
mainDir <- paste(rootDir, "/Data/Output", sep = "")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

# pdf file
cancertype <- subtype
outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
dat <- surdata[which(surdata$cancertype == cancertype),]

# Get subtype Basal
type <- "Basal"
coxdata <- as.data.frame(dat)
patients <- PAN.subtypes[which(PAN.subtypes$Subtype_mRNA == type),1]

coxdata <- coxdata[which(rownames(coxdata) %in% patients),]
coxdata[1:3,1:4]
nrow(coxdata)
# nrow(coxdata)
# cancertype OS OS.time  AANAT
# TCGA.A1.A0SK.01       BRCA  1     967 1.5850
# TCGA.A1.A0SO.01       BRCA  0     852 2.3219
# TCGA.A1.A0SP.01       BRCA  0     584 1.0000
# [1] 191

# Critical copper genes
criticalCopperGenes <- y
upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,"BRCA"] == 'UP'),]$gene
downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,"BRCA"] == 'DOWN'),]$gene
geneList <- criticalCopperGenes[which(criticalCopperGenes[,"BRCA"] %in% c("UP", "DOWN")),]$gene
fup <- paste(outDir, "/", cancertype, "_surUp_triple.pdf", sep = "")
fdown <- paste(outDir, "/", cancertype, "_surDown_triple.pdf", sep = "")
# > upGeneList
# [1] "CDK1"  "CASP3" "GSK3B" "XIAP"  "ARF1"  "SORD"
# > downGeneList
# [1] "ALB"      "MAP1LC3A" "SNCA"     "JUN"      "FOXO1"    "AOC3"     "S100A12"
# > geneList
# [1] "CDK1"     "ALB"      "CASP3"    "MAP1LC3A" "SNCA"     "GSK3B"
# [7] "JUN"      "XIAP"     "FOXO1"    "ARF1"     "AOC3"     "SORD"
# [13] "S100A12"

if(length(upGeneList) > 0) {
  pUp <- surAnalysis(coxdata, upGeneList, numGroup=2, K=20, alpha=0.5, fup, w = 9, h = 6)  
} else {
  pUp <- NA
}
if(length(downGeneList) > 0) {
  pDown <- surAnalysis(coxdata, downGeneList, numGroup=2, K=20, alpha=0.5, fdown, w = 9, h = 6)
} else {
  pDown <- NA
}

#================================================================

#================================================================
# (20) Compare gene expression, copper genes, 18 & 12 genes, 23 cancer types, using count and cpm
# run on server
#================================================================

# Link the files
# Katana
# cancertype="ACC"
# ln -s /srv/scratch/z3538133/001pancancer/pan/data/$cancertype/gtex.$cancertype.raw.counts.RData /srv/scratch/z3538133/002NetworkAnalysis/Data/$cancertype/gtex.$cancertype.raw.counts.RData

# Get the necessary data
getHeatmapData=function(geneList) {
  subtypes <- PanCancerAtlas_subtypes()
  subtypes <- unique(subtypes$cancer.type)
  subtypes <- c(subtypes, "PAAD")
  subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
  subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
  n <- length(subtypes)
  Cancer_type <- c()
  Type <- c()
  Expression <- c()
  for (i in 1:n) {
    print(i)
    cancertype <- subtypes[i]
    outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
    load(paste(outDir, "/gtex.", cancertype, ".raw.counts.RData", sep = "")) # return gtex.PAN.raw.counts, including both tumor and normal
    # # transform
    # dat <- gtex.PAN.raw.counts
    # rownames(dat) <- dat$gene
    # dat <- dat[,-1]
    # dat <- dat + 1
    # dat <- log(dat, base = 2)
    # iGenes <- which(rownames(dat) %in% geneList)
    # r <- str_detect(colnames(dat), "TCGA")
    # nTumor <- length(which(r == TRUE))
    # nNormal <- ncol(dat) - nTumor
    
    # compute logCPM
    dat <- gtex.PAN.raw.counts
    rownames(dat) <- dat$gene
    dat <- dat[,-1]
    # dat <- dat + 1
    logCPM <- cpm(dat, prior.count=2, log=TRUE)
    # scale columns, scale genes
    logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
    dat <- logCPM
    iGenes <- which(rownames(dat) %in% geneList)
    r <- str_detect(colnames(dat), "TCGA")
    nTumor <- length(which(r == TRUE))
    nNormal <- ncol(dat) - nTumor
    
    # # Add * if being critical copper gene
    # outDir <- paste(rootDir, "/Data/Output", sep = "")
    # fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
    # selectedCriticalCopperGenes <- read.csv(fileName)
    # if(selectedCriticalCopperGenes[which(selectedCriticalCopperGenes$gene == gene),cancertype] == 1){
    #   cancertype <- paste("*", cancertype, sep = "")
    # }
    
    Cancer_type <- c(Cancer_type, rep(cancertype, nTumor + nNormal))
    Type <- c(Type, rep(c("Tumour", "Normal"), times = c(nTumor, nNormal)))
    x <- dat[iGenes,]
    x <- t(x)
    # order by gene names
    x <- x[, order(colnames(x), decreasing = FALSE)]
    Expression <- rbind(Expression, x)
    
  }
  data=data.frame(Cancer_type, Type ,  Expression)
  
  return(data)
}

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# > head(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# > nrow(uni_gene)
# [1] 35

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
# > geneList18
# [1] "CDK1"    "AP1S1"   "CASP3"   "TMPRSS6" "GSK3B"   "APP"     "COX17"
# [8] "XIAP"    "ARF1"    "GPC1"    "SORD"    "ATP7A"   "SP1"     "MT-CO1"
# [15] "SLC11A2" "ATP6AP1" "ADAM10"  "CP"
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]
# > geneList12
# [1] "MAP1LC3A" "SNCA"     "MAPT"     "JUN"      "CYP1A1"   "AOC3"
# [7] "PRNP"     "DBH"      "S100A12"  "AQP1"     "MT1X"     "XAF1"

geneList18_12 <- c(geneList18, geneList12)
# > geneList18_12
# [1] "CDK1"     "AP1S1"    "CASP3"    "TMPRSS6"  "GSK3B"    "APP"
# [7] "COX17"    "XIAP"     "ARF1"     "GPC1"     "SORD"     "ATP7A"
# [13] "SP1"      "MT-CO1"   "SLC11A2"  "ATP6AP1"  "ADAM10"   "CP"
# [19] "MAP1LC3A" "SNCA"     "MAPT"     "JUN"      "CYP1A1"   "AOC3"
# [25] "PRNP"     "DBH"      "S100A12"  "AQP1"     "MT1X"     "XAF1"

data <- getHeatmapData(geneList18_12)
# Save a single object to a file
saveRDS(data, paste(outDir, "/dataForHeatmap.rds", sep = ""))

#-----------------------------------------
# For block of normal and tumour
# Restore it under a different name
my_data <- readRDS(paste(outDir, "/dataForHeatmap.rds", sep = ""))

# format the data
dat <- my_data
dat <- dat[,-c(1,2)] # remove columns of cancer type and type of normal or tumour
dat <- t(dat)
colnames(dat) <- paste(my_data[,1], "_", my_data[,2], sep = "")

# get mean value
r <- matrix(0, nrow = 30, ncol = 46)
rownames(r) <- rownames(dat)
colnames(r) <- c(1:46)
subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
n <- length(subtypes)
for (i in 1:n) {
  cancertype <- subtypes[i]
  tumour <- paste(cancertype, "_Tumour", sep = "")
  normal <-  paste(cancertype, "_Normal", sep = "")
  colnames(r)[i + 0] <- normal
  colnames(r)[i + 23] <- tumour
  for (j in 1:30) {
    r[j, i + 0] <- mean(dat[j, which(colnames(dat) == normal)])
    r[j, i + 23] <- mean(dat[j, which(colnames(dat) == tumour)])  
  }
  
}

# Divide by up genes and down genes
# geneList18 <- gsub("[-]",".",geneList18)
# geneList12 <- gsub("[-]",".",geneList12)
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% geneList18),]
rr <- rbind(rr,r[which(row.names(r) %in% geneList12),])
write.csv(rr, paste(outDir, "/rr.csv", sep = ""), row.names=TRUE)

# draw the chart
f <- paste(outDir, "/heatmap18_12.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, 18), rep(2,12))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, 23), rep(2,23))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#-----------------------------------------

#-----------------------------------------
# For pair of normal and tumour
# Restore it under a different name
my_data <- readRDS(paste(outDir, "/dataForHeatmap.rds", sep = ""))

# format the data
dat <- my_data
dat <- dat[,-c(1,2)] # remove columns of cancer type and type of normal or tumour
dat <- t(dat)
colnames(dat) <- paste(my_data[,1], "_", my_data[,2], sep = "")

# get mean value
r <- matrix(0, nrow = 30, ncol = 46)
rownames(r) <- rownames(dat)
colnames(r) <- c(1:46)
subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
n <- length(subtypes)
for (i in 1:n) {
  cancertype <- subtypes[i]
  tumour <- paste(cancertype, "_Tumour", sep = "")
  normal <-  paste(cancertype, "_Normal", sep = "")
  colnames(r)[(i-1)*2 + 1] <- normal
  colnames(r)[(i-1)*2 + 2] <- tumour
  for (j in 1:30) {
    r[j, (i-1)*2 + 1] <- mean(dat[j, which(colnames(dat) == normal)])
    r[j, (i-1)*2 + 2] <- mean(dat[j, which(colnames(dat) == tumour)])  
  }
  
}

# Divide by up genes and down genes
# geneList18 <- gsub("[-]",".",geneList18)
# geneList12 <- gsub("[-]",".",geneList12)
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% geneList18),]
rr <- rbind(rr,r[which(row.names(r) %in% geneList12),])
write.csv(rr, paste(outDir, "/rrForPair.csv", sep = ""), row.names=TRUE)

# draw the chart
f <- paste(outDir, "/heatmap18_12_pair.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, 18), rep(2,12))
rowSide <- brewer.pal(3, "Set1")[my_group]
# my_group2 <- c(rep(1, 23), rep(2,23))
# colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
# heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#-----------------------------------------

#================================================================

#================================================================
# (21) Plot gene networks
# run on server
#================================================================
subDir <- "network"
mainDir <- paste(rootDir, "/Data/Output", sep = "")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
n <- length(subtypes)

for (i in 1:n) {
  # i <- 11 # for LGG
  
  cancertype <- subtypes[i]
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  
  # Network
  network <- read.csv(paste(outDir, "/cancer_network.csv", sep = ""))
  
  # Critical copper genes
  load(paste(outDir, '/coxdata.Rdata', sep = "")) # return coxdata
  dat <- coxdata[,3:ncol(coxdata)]
  critical_copper_genes <- colnames(dat)
  
  # only get edges in which those genes regulate other genes
  refined_network1 <- network[which(network$cause %in% critical_copper_genes),]
  nodes <- c(refined_network1$cause, refined_network1$effect)
  nodes <- unique(nodes)
  n <- length(nodes)
  mat <- matrix(0, n, n)
  colnames(mat) <- nodes
  rownames(mat) <- nodes
  for (h in 1:nrow(refined_network1)) {
    mat[refined_network1$cause[h], refined_network1$effect[h]] <- 1
  }
  
  # plot
  net <- network(mat, directed = TRUE)
  net %v% "Gene" = ifelse(nodes[1:n] %in% critical_copper_genes, "Critical_copper_gene", "Ordinary_gene")
  net %v% "color" = ifelse(net %v% "Gene" == "Critical_copper_gene", "#F5AF58", "white")
  net %v% "size" = ifelse(net %v% "Gene"== "Critical_copper_gene", 10, 1)
  f <- paste(rootDir, "/Data/Output/network/", cancertype, "_gene_network.pdf", sep = "")
  pdf(file = f, width = 15, height =  10)
  # ggnet2(net,
  #        arrow.size = 3, arrow.gap = 0.02,
  #        # color = "Gene", palette = c("Ordinary_gene" = "#64BBE2", "Critical_copper_gene" = "#F5AF58"), 
  #        color = "Gene", palette = c("Ordinary_gene" = "white", "Critical_copper_gene" = "#F5AF58"), 
  #        size = "Gene", size.palette = c("Ordinary_gene" = 1, "Critical_copper_gene" = 10), 
  #        edge.size = 0.5, edge.color = "black",
  #        label = TRUE, label.size = 2.5)
  #        # label = critical_copper_genes, label.color = "black", label.size = 5, label.alpha = 0.75)
  print(ggnet2(net,
         arrow.size = 4, arrow.gap = 0.02,
         color = "color", 
         node.size = 12,
         edge.size = 0.5, edge.color = "black",
         label = TRUE, label.size = 2.5))
  dev.off()
}

#================================================================

#================================================================
# (22) Compare gene expression, copper genes, 18 & 12 genes, 23 cancer types, using tpm, unit log2(tpm+0.001)
# run on server
#================================================================

#----------------------------------------------------
# (22.1) For TCGA data
#----------------------------------------------------
# gzip -d tcga_RSEM_gene_tpm.gz # this will remove gz file
# gzip -d gtex_RSEM_gene_tpm.gz # this will remove gz file
# awk 'NR <= 5 {print $1, $2, $3, $4, $5, $6}' tcga_RSEM_gene_tpm
# (base) [z3538133@katana3 Data]$ awk 'NR <= 5 {print $1, $2, $3, $4, $5, $6}' tcga_RSEM_gene_tpm
# sample TCGA-19-1787-01 TCGA-S9-A7J2-01 TCGA-G3-A3CH-11 TCGA-EK-A2RE-01 TCGA-44-6778-01
# ENSG00000242268.2 -9.9658 0.2998 -9.9658 -9.9658 -9.9658
# ENSG00000259041.1 -9.9658 -9.9658 -9.9658 -9.9658 -9.9658
# ENSG00000270112.3 -3.8160 -3.0469 -9.9658 -9.9658 -5.5735
# ENSG00000167578.16 5.2998 4.8881 3.5572 4.2563 5.3162
PAN.file = paste(rootDir, "/Data/tcga_RSEM_gene_tpm.gz", sep ="")
PAN.log2.tmp = read.table(PAN.file, header = TRUE, stringsAsFactors = FALSE, row.names = "sample")
PAN.log2.tmp[1:5,1:6]
nrow(PAN.log2.tmp)
ncol(PAN.log2.tmp)
# > PAN.log2.tmp[1:5,1:6]
# TCGA.19.1787.01 TCGA.S9.A7J2.01 TCGA.G3.A3CH.11
# ENSG00000242268.2          -9.9658          0.2998         -9.9658
# ENSG00000259041.1          -9.9658         -9.9658         -9.9658
# ENSG00000270112.3          -3.8160         -3.0469         -9.9658
# ENSG00000167578.16          5.2998          4.8881          3.5572
# ENSG00000278814.1          -9.9658         -9.9658         -9.9658
# TCGA.EK.A2RE.01 TCGA.44.6778.01 TCGA.F4.6854.01
# ENSG00000242268.2          -9.9658         -9.9658         -9.9658
# ENSG00000259041.1          -9.9658         -9.9658         -9.9658
# ENSG00000270112.3          -9.9658         -5.5735         -9.9658
# ENSG00000167578.16          4.2563          5.3162          4.5161
# ENSG00000278814.1          -9.9658         -9.9658         -9.9658
# > nrow(PAN.log2.tmp)
# [1] 60498
# > ncol(PAN.log2.tmp)
# [1] 10535

# Create a subdataframe of data that only contains patient IDs from the PAN.path.subtypes dataframe
# subtypes <- PanCancerAtlas_subtypes()
# unique(subtypes$cancer.type)
# > unique(subtypes$cancer.type)
#  [1] "ACC"  "AML"  "BLCA" "BRCA" "LGG"  "GBM"  "ESCA" "COAD" "STAD" "READ" "HNSC" "KICH" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "OVCA" "PCPG" "PRAD"
# [21] "SKCM" "THCA" "UCEC" "UCS" 
# > nrow(subtypes)
# [1] 7734
subtypes <- c("ACC", "AML", "BLCA", "BRCA", "LGG", "GBM", "ESCA", "COAD", "STAD", "READ", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "OVCA", "PCPG", "PRAD", "SKCM", "THCA", "UCEC", "UCS")
all <- PanCancerAtlas_subtypes()
PAN.path.subtypes <- all[which(all$cancer.type %in% subtypes),]
PAN.subtypes <- PAN.path.subtypes
i <- sapply(PAN.subtypes, is.factor)
PAN.subtypes[i] <- lapply(PAN.subtypes[i], as.character)

# make list of column IDs
tumour.pt.ID <- as.list(PAN.subtypes$pan.samplesID) # makes list with 7734 elements

# replace - with . to match the colnames in the tpm dataframe
tmp1 <- str_replace_all(tumour.pt.ID, "-", ".")

# subset patient IDs that were only in the subtypes dataframe
# PAN.raw.counts <- PAN.raw.counts[,grep(paste(tmp1, collapse = "|"), x=names(PAN.raw.counts))] # Out of memory error
header <- names(PAN.log2.tmp)
header <- substr(header, 1, 12)
i <- which(header %in% tmp1)
PAN.log2.tmp <- PAN.log2.tmp[,i]
ncol(PAN.log2.tmp)
PAN.log2.tmp[1:5,1:6]
# > ncol(PAN.log2.tmp)
# [1] 3554
# > PAN.log2.tmp[1:5,1:6]
# TCGA.19.1787.01 TCGA.S9.A7J2.01 TCGA.F4.6854.01
# ENSG00000242268.2          -9.9658          0.2998         -9.9658
# ENSG00000259041.1          -9.9658         -9.9658         -9.9658
# ENSG00000270112.3          -3.8160         -3.0469         -9.9658
# ENSG00000167578.16          5.2998          4.8881          4.5161
# ENSG00000278814.1          -9.9658         -9.9658         -9.9658
# TCGA.AB.2863.03 TCGA.05.4420.01 TCGA.R6.A8WC.01
# ENSG00000242268.2          -9.9658         -9.9658         -9.9658
# ENSG00000259041.1          -9.9658         -9.9658         -9.9658
# ENSG00000270112.3          -3.4580         -9.9658         -6.5064
# ENSG00000167578.16          3.6242          4.2048          3.8581
# ENSG00000278814.1          -9.9658         -9.9658         -9.9658

# Convert ensembl IDs to gene names
genecode.file <- paste(rootDir, "/Data/probeMap_gencode.v23.annotation.gene.probemap", sep = "")
genecode <- read.table(genecode.file, header = TRUE, sep = "\t")
genecode.tmp <- genecode[, 1:2]

# convert rownames in dataframe into a column
PAN.log2.tmp$Ensembl_ID <- rownames(PAN.log2.tmp)
PAN.log2.tmp <- PAN.log2.tmp %>% dplyr::select(Ensembl_ID, everything())

# match ensembl ID of genecode dataframe to PAN.log2.tmp dataframe
PAN.log2.tmp <- merge(PAN.log2.tmp, genecode.tmp, by.x = "Ensembl_ID", by.y = "id")
PAN.log2.tmp <- PAN.log2.tmp %>% dplyr::select(gene, everything())
PAN.log2.tmp <- PAN.log2.tmp[ ,-2] # remove Ensembl_ID
PAN.log2.tmp.unique <- aggregate(PAN.log2.tmp, by=list(PAN.log2.tmp$gene), FUN=mean)
PAN.log2.tmp.unique[1:5,1:6]
# > PAN.log2.tmp.unique[1:5,1:6]
# Group.1 gene TCGA.19.1787.01 TCGA.S9.A7J2.01 TCGA.F4.6854.01
# 1 5_8S_rRNA   NA         -9.9658         -9.9658         -9.9658
# 2   5S_rRNA   NA         -9.9658         -9.9658         -9.9658
# 3       7SK   NA         -9.9658         -9.9658         -9.9658
# 4      A1BG   NA          2.0810          4.4842          1.2085
# 5  A1BG-AS1   NA          0.9493          2.0707         -0.8599
# TCGA.AB.2863.03
# 1       -9.965800
# 2       -9.081019
# 3       -6.794567
# 4        2.060400
# 5        2.696200
PAN.log2.tmp.unique$gene <- PAN.log2.tmp.unique$Group.1
PAN.log2.tmp.unique <- PAN.log2.tmp.unique[ ,-1]
# > PAN.log2.tmp.unique[1:5,1:6]
# gene TCGA.19.1787.01 TCGA.S9.A7J2.01 TCGA.F4.6854.01 TCGA.AB.2863.03
# 1 5_8S_rRNA         -9.9658         -9.9658         -9.9658       -9.965800
# 2   5S_rRNA         -9.9658         -9.9658         -9.9658       -9.081019
# 3       7SK         -9.9658         -9.9658         -9.9658       -6.794567
# 4      A1BG          2.0810          4.4842          1.2085        2.060400
# 5  A1BG-AS1          0.9493          2.0707         -0.8599        2.696200
# TCGA.05.4420.01
# 1         -9.9658
# 2         -9.9658
# 3         -9.9658
# 4          5.9298
# 5          1.7912

# save data
# write.csv(PAN.raw.counts.unique, paste(rootDir, "/data/output/",'PAN.raw.counts.unique.csv', sep = ""), row.names = F)
save(PAN.log2.tmp.unique, file=paste(rootDir, "/Data/Output/",'PAN.log2.tmp.unique.RData', sep = ""))  

#----------------------------------------------------

# #----------------------------------------------------
# # (22.2) For GTEX data
# #----------------------------------------------------
# 
# GTEX.file <- paste(rootDir, "/Data/gtex_RSEM_gene_tpm.gz", sep = "")
# gtex.log2.tmp <- read.table(GTEX.file, header = TRUE)
# # > gtex.log2.tmp[1:5,1:6]
# # sample GTEX.S4Q7.0003.SM.3NM8M GTEX.QV31.1626.SM.2S1QC
# # 1  ENSG00000242268.2                 -3.4580                 -9.9658
# # 2  ENSG00000259041.1                 -9.9658                 -9.9658
# # 3  ENSG00000270112.3                 -3.6259                 -2.1779
# # 4 ENSG00000167578.16                  4.5988                  4.6294
# # 5  ENSG00000278814.1                 -9.9658                 -9.9658
# # GTEX.13QIC.0011.R1a.SM.5O9CJ GTEX.ZPCL.0126.SM.4WWC8 GTEX.S33H.1226.SM.4AD69
# # 1                      -9.9658                 -9.9658                 -2.7274
# # 2                      -9.9658                 -9.9658                 -9.9658
# # 3                      -1.8314                 -9.9658                 -9.9658
# # 4                       6.4989                  5.5358                  3.7269
# # 5                      -9.9658                 -9.9658                 -9.9658
# # > nrow(gtex.log2.tmp)
# # [1] 60498
# # > ncol(gtex.log2.tmp)
# # [1] 7863
# GTEX.pheno.file <- paste(rootDir, "/Data/GTEX_phenotype.gz", sep = "")
# gtex.pheno <- read.table(GTEX.pheno.file, header = TRUE, sep = "\t")
# # > gtex.pheno[1:5,1:6]
# # Sample body_site_detail..SMTSD. X_primary_site X_gender
# # 1 GTEX-1117F-0226-SM-5GZZ7   Adipose - Subcutaneous Adipose Tissue   female
# # 2 GTEX-1117F-0426-SM-5EGHI        Muscle - Skeletal         Muscle   female
# # 3 GTEX-1117F-0526-SM-5EGHJ          Artery - Tibial   Blood Vessel   female
# # 4 GTEX-1117F-0626-SM-5N9CS        Artery - Coronary   Blood Vessel   female
# # 5 GTEX-1117F-0726-SM-5GIEN Heart - Atrial Appendage          Heart   female
# # X_patient X_cohort
# # 1 GTEX-1117F     GTEX
# # 2 GTEX-1117F     GTEX
# # 3 GTEX-1117F     GTEX
# # 4 GTEX-1117F     GTEX
# # 5 GTEX-1117F     GTEX
# # > nrow(gtex.pheno)
# # [1] 9783
# # > ncol(gtex.pheno)
# # [1] 6
# 
# # Create subset dataframe to contain only disease samples
# # first filter out other rows not in disease
# # > unique(gtex.pheno$X_primary_site)
# #  [1] "Adipose Tissue"  "Muscle"          "Blood Vessel"    "Heart"           "Ovary"           "Uterus"          "Vagina"          "Breast"          "Skin"           
# # [10] "Salivary Gland"  "Brain"           "Adrenal Gland"   "Thyroid"         "Lung"            "Spleen"          "Pancreas"        "Esophagus"       "Stomach"        
# # [19] "Colon"           "Small Intestine" "Prostate"        "Testis"          "Nerve"           "Pituitary"       "Blood"           "Liver"           "Kidney"         
# # [28] "Fallopian Tube"  "Bladder"         "Cervix Uteri"    "<not provided>"  "Bone Marrow"
# gtex <- gtex.pheno[gtex.pheno$X_primary_site != "<not provided>", ]
# # > nrow(gtex)
# # [1] 9778
# 
# # make list of column IDs
# gtex.sample.ID <- as.list(gtex$Sample)
# gtex.sample.ID <- append(gtex.sample.ID, "sample")
# 
# # replace - with . to match the colnames in the dataframe
# tmp2 <- str_replace_all(gtex.sample.ID, "-", ".")
# 
# # subset patient IDs that were only in the subtypes dataframe
# # pan.log2.counts <- gtex.log2.counts[,grep(paste(tmp2, collapse = "|"), x=names(gtex.log2.counts))]
# header <- names(gtex.log2.tmp)
# i <- which(header %in% tmp2)
# # > head(i)
# # [1] 1 2 3 4 5 6
# # > length(i)
# # [1] 7859
# pan.log2.tmp <- gtex.log2.tmp[,i]
# # > pan.log2.tmp[1:5,1:6]
# # sample GTEX.S4Q7.0003.SM.3NM8M GTEX.QV31.1626.SM.2S1QC
# # 1  ENSG00000242268.2                 -3.4580                 -9.9658
# # 2  ENSG00000259041.1                 -9.9658                 -9.9658
# # 3  ENSG00000270112.3                 -3.6259                 -2.1779
# # 4 ENSG00000167578.16                  4.5988                  4.6294
# # 5  ENSG00000278814.1                 -9.9658                 -9.9658
# # GTEX.13QIC.0011.R1a.SM.5O9CJ GTEX.ZPCL.0126.SM.4WWC8 GTEX.S33H.1226.SM.4AD69
# # 1                      -9.9658                 -9.9658                 -2.7274
# # 2                      -9.9658                 -9.9658                 -9.9658
# # 3                      -1.8314                 -9.9658                 -9.9658
# # 4                       6.4989                  5.5358                  3.7269
# # 5                      -9.9658                 -9.9658                 -9.9658
# # > nrow(pan.log2.tmp)
# # [1] 60498
# # > ncol(pan.log2.tmp)
# # [1] 7859
# 
# # assign gene column as rowname
# rownames(pan.log2.tmp) <- pan.log2.tmp[,"sample"]
# pan.log2.tmp[ ,"sample"] <- NULL
# # > pan.log2.tmp[1:5,1:4]
# # GTEX.S4Q7.0003.SM.3NM8M GTEX.QV31.1626.SM.2S1QC
# # ENSG00000242268.2                  -3.4580                 -9.9658
# # ENSG00000259041.1                  -9.9658                 -9.9658
# # ENSG00000270112.3                  -3.6259                 -2.1779
# # ENSG00000167578.16                  4.5988                  4.6294
# # ENSG00000278814.1                  -9.9658                 -9.9658
# # GTEX.13QIC.0011.R1a.SM.5O9CJ GTEX.ZPCL.0126.SM.4WWC8
# # ENSG00000242268.2                       -9.9658                 -9.9658
# # ENSG00000259041.1                       -9.9658                 -9.9658
# # ENSG00000270112.3                       -1.8314                 -9.9658
# # ENSG00000167578.16                       6.4989                  5.5358
# # ENSG00000278814.1                       -9.9658                 -9.9658
# 
# # convert ensembl IDs to gene names
# gtex.genecode.file <- paste(rootDir, "/Data/probeMap_gencode.v23.annotation.gene.probemap", sep = "")
# gtex.genecode <- read.table(gtex.genecode.file, header = TRUE, sep = "\t")
# # > head(gtex.genecode)
# # id         gene chrom chromStart chromEnd strand
# # 1 ENSG00000223972.5      DDX11L1  chr1      11869    14409      +
# #   2 ENSG00000227232.5       WASH7P  chr1      14404    29570      -
# #   3 ENSG00000278267.1    MIR6859-1  chr1      17369    17436      -
# #   4 ENSG00000243485.3 RP11-34P13.3  chr1      29554    31109      +
# #   5 ENSG00000274890.1    MIR1302-2  chr1      30366    30503      +
# #   6 ENSG00000237613.2      FAM138A  chr1      34554    36081      -
# gtex.genecode.tmp <- gtex.genecode[, 1:2]
# 
# # assign rownames as a column so you can use the merge() function
# pan.log2.tmp <- tibble::rownames_to_column(pan.log2.tmp, "Ensembl_ID") # Apply rownames_to_column
# 
# # match ensembl ID of genecode dataframe to log2.tmp dataframe
# pan.log2.tmp <- merge(pan.log2.tmp, gtex.genecode.tmp, by.x = "Ensembl_ID", by.y = "id")
# pan.log2.tmp <- pan.log2.tmp %>% dplyr::select(gene, everything())
# pan.log2.tmp <- pan.log2.tmp[ ,-2]
# pan.log2.tmp.unique <- aggregate(pan.log2.tmp, by=list(pan.log2.tmp$gene), FUN=mean)
# pan.log2.tmp.unique$gene <-pan.log2.tmp.unique$Group.1
# pan.log2.tmp.unique <- pan.log2.tmp.unique[ ,-1]
# 
# # save data
# write.csv(pan.log2.tmp.unique,paste(rootDir, "/Data/Output/",'pan.log2.tmp.unique.csv', sep = ""), row.names = F)
# save(pan.log2.tmp.unique, file=paste(rootDir, "/Data/Output/",'pan.log2.tmp.unique.RData', sep = ""))
# 
# #----------------------------------------------------

#================================================================

#================================================================
# (23) Heat map for Survival analysis, using 18 genes and 12 genes for each cancer type
# Run in server
#================================================================

# Get the necessary data for the heat map for each cancer type
# rootDir: Root directory, the data will be in rootDir/Data/...
# geneList: List of genes we need the data
# cancertype: The cancer for the analysis
getHeatmapDataForCancerType=function(rootDir, geneList, cancertype) {
  
  Cancer_type <- c()
  Type <- c()
  Expression <- c()
  
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(paste(outDir, "/gtex.", cancertype, ".raw.counts.RData", sep = "")) # return gtex.PAN.raw.counts, including both tumor and normal
  # # transform
  # dat <- gtex.PAN.raw.counts
  # rownames(dat) <- dat$gene
  # dat <- dat[,-1]
  # dat <- dat + 1
  # dat <- log(dat, base = 2)
  # iGenes <- which(rownames(dat) %in% geneList)
  # r <- str_detect(colnames(dat), "TCGA")
  # nTumor <- length(which(r == TRUE))
  # nNormal <- ncol(dat) - nTumor
  
  # compute logCPM
  dat <- gtex.PAN.raw.counts
  rownames(dat) <- dat$gene
  dat <- dat[,-1]
  # dat <- dat + 1
  logCPM <- cpm(dat, prior.count=2, log=TRUE)
  logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
  dat <- logCPM
  iGenes <- which(rownames(dat) %in% geneList)
  r <- str_detect(colnames(dat), "TCGA")
  nTumor <- length(which(r == TRUE))
  nNormal <- ncol(dat) - nTumor
  
  # # Add * if being critical copper gene
  # outDir <- paste(rootDir, "/Data/Output", sep = "")
  # fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
  # selectedCriticalCopperGenes <- read.csv(fileName)
  # if(selectedCriticalCopperGenes[which(selectedCriticalCopperGenes$gene == gene),cancertype] == 1){
  #   cancertype <- paste("*", cancertype, sep = "")
  # }
  
  Cancer_type <- c(Cancer_type, rep(cancertype, nTumor + nNormal))
  Type <- c(Type, rep(c("Tumour", "Normal"), times = c(nTumor, nNormal)))
  x <- dat[iGenes,]
  x <- t(x)
  # order by gene names
  x <- x[, order(colnames(x), decreasing = FALSE)]
  Expression <- rbind(Expression, x)
  
  data=data.frame(Cancer_type, Type ,  Expression)
  
  return(data)
}

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# > head(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# > nrow(uni_gene)
# [1] 35

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
# > geneList18
# [1] "CDK1"    "AP1S1"   "CASP3"   "TMPRSS6" "GSK3B"   "APP"     "COX17"
# [8] "XIAP"    "ARF1"    "GPC1"    "SORD"    "ATP7A"   "SP1"     "MT-CO1"
# [15] "SLC11A2" "ATP6AP1" "ADAM10"  "CP"
# f18 <- paste(outDir, "/sur18.pdf", sep = "")
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]
# > geneList12
# [1] "MAP1LC3A" "SNCA"     "MAPT"     "JUN"      "CYP1A1"   "AOC3"
# [7] "PRNP"     "DBH"      "S100A12"  "AQP1"     "MT1X"     "XAF1"
# f12 <- paste(outDir, "/sur12.pdf", sep = "")
# K	Number of nearest neighbors, default is 20
# surAnalysis(surdata, geneList18, numGroup=2, K=20, alpha=0.5, f18, w = 9, h = 6)
# surAnalysis(surdata, geneList12, numGroup=2, K=20, alpha=0.5, f12, w = 9, h = 6)
print(paste(geneList18, collapse = ", "))
print(paste(geneList12, collapse = ", "))
# > print(paste(geneList18, collapse = ", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, ARF1, GPC1, SORD, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"
# > print(paste(geneList12, collapse = ", "))
# [1] "MAP1LC3A, SNCA, MAPT, JUN, CYP1A1, AOC3, PRNP, DBH, S100A12, AQP1, MT1X, XAF1"

geneList18_12 <- c(geneList18, geneList12)
# > geneList18_12
# [1] "CDK1"     "AP1S1"    "CASP3"    "TMPRSS6"  "GSK3B"    "APP"
# [7] "COX17"    "XIAP"     "ARF1"     "GPC1"     "SORD"     "ATP7A"
# [13] "SP1"      "MT-CO1"   "SLC11A2"  "ATP6AP1"  "ADAM10"   "CP"
# [19] "MAP1LC3A" "SNCA"     "MAPT"     "JUN"      "CYP1A1"   "AOC3"
# [25] "PRNP"     "DBH"      "S100A12"  "AQP1"     "MT1X"     "XAF1"

data <- getHeatmapDataForCancerType(rootDir, geneList18_12, "LGG")
# Save a single object to a file
cancertype <- "LGG"
saveRDS(data, paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))

#---------------------------------------
# For 18 & 12 genes, for normal and tumour
# Restore it under a different name
my_data <- readRDS(paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))
# > my_data[1:4,1:5]
# Cancer_type   Type    ADAM10        AOC3      AP1S1
# TCGA.S9.A7J2.01         LGG Tumour 0.4730344 -0.22753447  0.1251829
# TCGA.DU.7302.01         LGG Tumour 0.9706330 -1.19686021  0.8742207
# TCGA.P5.A5EX.01         LGG Tumour 0.3353221  0.79421141 -0.9714453
# TCGA.HT.8563.01         LGG Tumour 1.1853219  0.09386498 -0.8300805
# > nrow(my_data)
# [1] 453
# > ncol(my_data)
# [1] 32
# my_data <- my_data[which(my_data$Type == "Tumour"),]
my_data$patient_id <- rownames(my_data)
my_data <- my_data %>% dplyr::select(patient_id, everything())
# > my_data[1:4,1:5]
# patient_id Cancer_type   Type    ADAM10        AOC3
# TCGA.S9.A7J2.01 TCGA.S9.A7J2.01         LGG Tumour 0.4730344 -0.22753447
# TCGA.DU.7302.01 TCGA.DU.7302.01         LGG Tumour 0.9706330 -1.19686021
# TCGA.P5.A5EX.01 TCGA.P5.A5EX.01         LGG Tumour 0.3353221  0.79421141
# TCGA.HT.8563.01 TCGA.HT.8563.01         LGG Tumour 1.1853219  0.09386498
# > nrow(my_data)
# [1] 453

outDir <- paste(rootDir, "/Data/Output/surCancerType1812", sep = "")
# load(paste(outDir, "/", cancertype, "_sur18.Rdata", sep = "")) # return result1
# 
# patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
# colnames(patients) <- c("patient_id", "group")
# patients[,1] <- as.character(colnames(result1$distanceMatrix))
# patients[,2] <- result1$group
# # > head(patients)
# # patient_id        group
# # [1,] "TCGA.CS.4938.01" "2"
# # [2,] "TCGA.CS.4941.01" "2"
# # [3,] "TCGA.CS.4942.01" "2"
# # [4,] "TCGA.CS.4943.01" "2"
# # [5,] "TCGA.CS.4944.01" "2"
# # [6,] "TCGA.CS.5390.01" "2"

# new_dat <- merge(my_data, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
# new_dat <- new_dat %>% dplyr::select(group, everything())
# new_dat <- new_dat[order(new_dat$group),]
# # > new_dat[1:4,1:5]
# # group      patient_id Cancer_type   Type    ADAM10
# # 7      1 TCGA.CS.5393.01         LGG Tumour 0.8411652
# # 8      1 TCGA.CS.5394.01         LGG Tumour 1.0374810
# # 9      1 TCGA.CS.5395.01         LGG Tumour 0.4455639
# # 15     1 TCGA.CS.6667.01         LGG Tumour 0.1035104
# # > nrow(new_dat)
# # [1] 348
# # > ncol(new_dat)
# # [1] 34

# format the data
# dat <- new_dat
dat <- my_data
dat <- dat[,-c(1,2,3)]
dat <- t(dat)
colnames(dat) <- paste(my_data[,1], "_", my_data[,3], sep = "")

# For 18 genes & 12 genes
r <- dat
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% geneList18),]
rr <- rbind(rr, r[which(row.names(r) %in% geneList12),])
nTumour <- nrow(my_data[which(my_data$Type == "Tumour"),])
nNormal <- nrow(my_data[which(my_data$Type == "Normal"),])
tmp <- rr
rr <- tmp[,(nTumour+1):(nTumour + nNormal)]
rr <- cbind(rr, tmp[,1:nTumour])
# draw
f <- paste(outDir, "/", cancertype, "_heatmap18_12.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, 18), rep(2,12))
# my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, nNormal), rep(2,nTumour))
# my_group2 <- c(rep(1, nrow(new_dat[which(new_dat$group == 1),])),
#                rep(2, nrow(new_dat[which(new_dat$group == 2),])))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#---------------------------------------

#---------------------------------------
# For 18 genes, for 2 groups

# Restore it under a different name
my_data <- readRDS(paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))
# > my_data[1:4,1:5]
# Cancer_type   Type    ADAM10        AOC3      AP1S1
# TCGA.S9.A7J2.01         LGG Tumour 0.4730344 -0.22753447  0.1251829
# TCGA.DU.7302.01         LGG Tumour 0.9706330 -1.19686021  0.8742207
# TCGA.P5.A5EX.01         LGG Tumour 0.3353221  0.79421141 -0.9714453
# TCGA.HT.8563.01         LGG Tumour 1.1853219  0.09386498 -0.8300805
# > nrow(my_data)
# [1] 453
# > ncol(my_data)
# [1] 32
my_data <- my_data[which(my_data$Type == "Tumour"),]
my_data$patient_id <- rownames(my_data)
my_data <- my_data %>% dplyr::select(patient_id, everything())
# > my_data[1:4,1:5]
# patient_id Cancer_type   Type    ADAM10        AOC3
# TCGA.S9.A7J2.01 TCGA.S9.A7J2.01         LGG Tumour 0.4730344 -0.22753447
# TCGA.DU.7302.01 TCGA.DU.7302.01         LGG Tumour 0.9706330 -1.19686021
# TCGA.P5.A5EX.01 TCGA.P5.A5EX.01         LGG Tumour 0.3353221  0.79421141
# TCGA.HT.8563.01 TCGA.HT.8563.01         LGG Tumour 1.1853219  0.09386498
# > nrow(my_data)
# [1] 348

outDir <- paste(rootDir, "/Data/Output/surCancerType1812", sep = "")
load(paste(outDir, "/", cancertype, "_sur18.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "2"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "2"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "2"
# [6,] "TCGA.CS.5390.01" "2"

new_dat <- merge(my_data, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
new_dat <- new_dat %>% dplyr::select(group, everything())
new_dat <- new_dat[order(new_dat$group),]
# > new_dat[1:4,1:5]
# group      patient_id Cancer_type   Type    ADAM10
# 7      1 TCGA.CS.5393.01         LGG Tumour 0.8411652
# 8      1 TCGA.CS.5394.01         LGG Tumour 1.0374810
# 9      1 TCGA.CS.5395.01         LGG Tumour 0.4455639
# 15     1 TCGA.CS.6667.01         LGG Tumour 0.1035104
# > nrow(new_dat)
# [1] 348
# > ncol(new_dat)
# [1] 34

# format the data
dat <- new_dat
# dat <- my_data
dat <- dat[,-c(1,2,3,4)]
colnames(dat) <- gsub("[.]","-",colnames(dat))
dat <- dat[, which(colnames(dat) %in% geneList18)]
dat <- t(dat)
colnames(dat) <- paste(new_dat[,2], "_", new_dat[,1], sep = "")
# > dat[1:5,1:4]
# TCGA.CS.5393.01_1 TCGA.CS.5394.01_1 TCGA.CS.5395.01_1 TCGA.CS.6667.01_1
# ADAM10          0.8411652         1.0374810        0.44556385        0.10351044
# AP1S1          -0.8580676         0.9990138       -0.08959111       -1.67823550
# APP            -0.2354364         0.5455580       -0.38693952       -0.08129716
# ARF1           -0.6763754        -0.9576393       -1.61662152       -0.20378669
# ATP6AP1        -0.8344859        -0.2209747       -1.69036999       -0.08097393
# > nrow(dat)
# [1] 18
# > ncol(dat)
# [1] 348

# For 18 genes
r <- dat
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% geneList18),]
#rr <- rbind(rr, r[which(row.names(r) %in% geneList12),])
nGroup1 <- nrow(patients[which(patients[,2] == 1),])
nGroup2 <- nrow(patients[which(patients[,2] == 2),])
# draw
f <- paste(outDir, "/", cancertype, "_heatmap18.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, 18))
# my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, nGroup1), rep(2,nGroup2))
# my_group2 <- c(rep(1, nrow(new_dat[which(new_dat$group == 1),])),
#                rep(2, nrow(new_dat[which(new_dat$group == 2),])))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

#---------------------------------------

#---------------------------------------
# For 12 genes, for 2 groups

# Restore it under a different name
my_data <- readRDS(paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))
# > my_data[1:4,1:5]
# Cancer_type   Type    ADAM10        AOC3      AP1S1
# TCGA.S9.A7J2.01         LGG Tumour 0.4730344 -0.22753447  0.1251829
# TCGA.DU.7302.01         LGG Tumour 0.9706330 -1.19686021  0.8742207
# TCGA.P5.A5EX.01         LGG Tumour 0.3353221  0.79421141 -0.9714453
# TCGA.HT.8563.01         LGG Tumour 1.1853219  0.09386498 -0.8300805
# > nrow(my_data)
# [1] 453
# > ncol(my_data)
# [1] 32
my_data <- my_data[which(my_data$Type == "Tumour"),]
my_data$patient_id <- rownames(my_data)
my_data <- my_data %>% dplyr::select(patient_id, everything())
# > my_data[1:4,1:5]
# patient_id Cancer_type   Type    ADAM10        AOC3
# TCGA.S9.A7J2.01 TCGA.S9.A7J2.01         LGG Tumour 0.4730344 -0.22753447
# TCGA.DU.7302.01 TCGA.DU.7302.01         LGG Tumour 0.9706330 -1.19686021
# TCGA.P5.A5EX.01 TCGA.P5.A5EX.01         LGG Tumour 0.3353221  0.79421141
# TCGA.HT.8563.01 TCGA.HT.8563.01         LGG Tumour 1.1853219  0.09386498
# > nrow(my_data)
# [1] 348

outDir <- paste(rootDir, "/Data/Output/surCancerType1812", sep = "")
load(paste(outDir, "/", cancertype, "_sur12.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "2"
# [2,] "TCGA.CS.4941.01" "1"
# [3,] "TCGA.CS.4942.01" "2"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "2"
# [6,] "TCGA.CS.5390.01" "2"

new_dat <- merge(my_data, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
new_dat <- new_dat %>% dplyr::select(group, everything())
new_dat <- new_dat[order(new_dat$group),]
# > new_dat[1:4,1:5]
# group      patient_id Cancer_type   Type    ADAM10
# 2      1 TCGA.CS.4941.01         LGG Tumour 0.3482559
# 9      1 TCGA.CS.5395.01         LGG Tumour 0.4455639
# 11     1 TCGA.CS.5397.01         LGG Tumour 0.5036000
# 12     1 TCGA.CS.6188.01         LGG Tumour 1.6559567
# > nrow(new_dat)
# [1] 348
# > ncol(new_dat)
# [1] 34

# format the data
dat <- new_dat
# dat <- my_data
dat <- dat[,-c(1,2,3,4)]
colnames(dat) <- gsub("[.]","-",colnames(dat))
dat <- dat[, which(colnames(dat) %in% geneList12)]
dat <- t(dat)
colnames(dat) <- paste(new_dat[,2], "_", new_dat[,1], sep = "")
# > dat[1:5,1:4]
# TCGA.CS.4941.01_1 TCGA.CS.5395.01_1 TCGA.CS.5397.01_1 TCGA.CS.6188.01_1
# AOC3          -1.3205445         2.0679927        -0.3313792        -0.4341555
# AQP1           2.2255065         2.3324790         1.6683529         1.5165818
# CYP1A1        -1.5145701        -1.1529441         0.3998401        -1.2926878
# DBH           -0.8451963         0.3102576         0.3214936        -0.4322281
# JUN            0.9109134        -0.1858019         0.2628002         0.4970608
# > nrow(dat)
# [1] 12
# > ncol(dat)
# [1] 348

# For 12 genes
r <- dat
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% geneList12),]
#rr <- rbind(rr, r[which(row.names(r) %in% geneList12),])
nGroup1 <- nrow(patients[which(patients[,2] == 1),])
nGroup2 <- nrow(patients[which(patients[,2] == 2),])
# draw
f <- paste(outDir, "/", cancertype, "_heatmap12.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, 12))
# my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, nGroup1), rep(2,nGroup2))
# my_group2 <- c(rep(1, nrow(new_dat[which(new_dat$group == 1),])),
#                rep(2, nrow(new_dat[which(new_dat$group == 2),])))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#---------------------------------------

#================================================================

#================================================================
# (24) Heat map for Survival analysis, using critical copper genes for each cancer type
# Run in server
#================================================================

# Get the necessary data for the heat map for each cancer type
# rootDir: Root directory, the data will be in rootDir/Data/...
# geneList: List of genes we need the data
# cancertype: The cancer for the analysis
getHeatmapDataForCancerType=function(rootDir, geneList, cancertype) {
  
  Cancer_type <- c()
  Type <- c()
  Expression <- c()
  
  outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
  load(paste(outDir, "/gtex.", cancertype, ".raw.counts.RData", sep = "")) # return gtex.PAN.raw.counts, including both tumor and normal
  # # transform
  # dat <- gtex.PAN.raw.counts
  # rownames(dat) <- dat$gene
  # dat <- dat[,-1]
  # dat <- dat + 1
  # dat <- log(dat, base = 2)
  # iGenes <- which(rownames(dat) %in% geneList)
  # r <- str_detect(colnames(dat), "TCGA")
  # nTumor <- length(which(r == TRUE))
  # nNormal <- ncol(dat) - nTumor
  
  # compute logCPM
  dat <- gtex.PAN.raw.counts
  rownames(dat) <- dat$gene
  dat <- dat[,-1]
  # dat <- dat + 1
  logCPM <- cpm(dat, prior.count=2, log=TRUE)
  logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
  dat <- logCPM
  iGenes <- which(rownames(dat) %in% geneList)
  r <- str_detect(colnames(dat), "TCGA")
  nTumor <- length(which(r == TRUE))
  nNormal <- ncol(dat) - nTumor
  
  # # Add * if being critical copper gene
  # outDir <- paste(rootDir, "/Data/Output", sep = "")
  # fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
  # selectedCriticalCopperGenes <- read.csv(fileName)
  # if(selectedCriticalCopperGenes[which(selectedCriticalCopperGenes$gene == gene),cancertype] == 1){
  #   cancertype <- paste("*", cancertype, sep = "")
  # }
  
  Cancer_type <- c(Cancer_type, rep(cancertype, nTumor + nNormal))
  Type <- c(Type, rep(c("Tumour", "Normal"), times = c(nTumor, nNormal)))
  x <- dat[iGenes,]
  x <- t(x)
  # order by gene names
  x <- x[, order(colnames(x), decreasing = FALSE)]
  Expression <- rbind(Expression, x)
  
  data=data.frame(Cancer_type, Type ,  Expression)
  
  return(data)
}

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
# > y[1:4,1:5]
# gene frequency ACC  AML BLCA
# 1  CDK1        18  UP   UP   UP
# 2   ALB        13   0 DOWN    0
# 3 AP1S1        11   0 DOWN    0
# 4 CASP3         9   0    0    0

subDir <- "surCancerTypeUseCopperGenes"
mainDir <- paste(rootDir, "/Data/Output", sep = "")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)

cancertype <- "LGG"

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")

# Critical copper genes
criticalCopperGenes <- y
upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'UP'),]$gene
downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'DOWN'),]$gene
geneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] %in% c("UP", "DOWN")),]$gene

data <- getHeatmapDataForCancerType(rootDir, geneList, "LGG")
# Save a single object to a file
saveRDS(data, paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))

#---------------------------------------
# For copper genes, for normal and tumour
# Restore it under a different name
my_data <- readRDS(paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))
# > my_data[1:4,1:5]
# Cancer_type   Type        ALB     ATP7A       CASP3
# TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519  1.10979667
# TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239 -0.67663103
# TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562  0.01558295
# TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199  0.79229205
# > nrow(my_data)
# [1] 453
# > ncol(my_data)
# [1] 15
# my_data <- my_data[which(my_data$Type == "Tumour"),]
my_data$patient_id <- rownames(my_data)
my_data <- my_data %>% dplyr::select(patient_id, everything())
# > my_data[1:4,1:5]
# patient_id Cancer_type   Type        ALB     ATP7A
# TCGA.S9.A7J2.01 TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519
# TCGA.DU.7302.01 TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239
# TCGA.P5.A5EX.01 TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562
# TCGA.HT.8563.01 TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199
# > nrow(my_data)
# [1] 453

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
# load(paste(outDir, "/", cancertype, "_sur18.Rdata", sep = "")) # return result1
# 
# patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
# colnames(patients) <- c("patient_id", "group")
# patients[,1] <- as.character(colnames(result1$distanceMatrix))
# patients[,2] <- result1$group
# # > head(patients)
# # patient_id        group
# # [1,] "TCGA.CS.4938.01" "2"
# # [2,] "TCGA.CS.4941.01" "2"
# # [3,] "TCGA.CS.4942.01" "2"
# # [4,] "TCGA.CS.4943.01" "2"
# # [5,] "TCGA.CS.4944.01" "2"
# # [6,] "TCGA.CS.5390.01" "2"

# new_dat <- merge(my_data, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
# new_dat <- new_dat %>% dplyr::select(group, everything())
# new_dat <- new_dat[order(new_dat$group),]
# # > new_dat[1:4,1:5]
# # group      patient_id Cancer_type   Type    ADAM10
# # 7      1 TCGA.CS.5393.01         LGG Tumour 0.8411652
# # 8      1 TCGA.CS.5394.01         LGG Tumour 1.0374810
# # 9      1 TCGA.CS.5395.01         LGG Tumour 0.4455639
# # 15     1 TCGA.CS.6667.01         LGG Tumour 0.1035104
# # > nrow(new_dat)
# # [1] 348
# # > ncol(new_dat)
# # [1] 34

# format the data
# dat <- new_dat
dat <- my_data
dat <- dat[,-c(1,2,3)]
dat <- t(dat)
colnames(dat) <- paste(my_data[,1], "_", my_data[,3], sep = "")

# For copper genes
r <- dat
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% upGeneList),]
rr <- rbind(rr, r[which(row.names(r) %in% downGeneList),])
nTumour <- nrow(my_data[which(my_data$Type == "Tumour"),])
nNormal <- nrow(my_data[which(my_data$Type == "Normal"),])
tmp <- rr
rr <- tmp[,(nTumour+1):(nTumour + nNormal)]
rr <- cbind(rr, tmp[,1:nTumour])
# draw
f <- paste(outDir, "/", cancertype, "_heatmapAll.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, length(upGeneList)), rep(2,length(downGeneList)))
# my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, nNormal), rep(2,nTumour))
# my_group2 <- c(rep(1, nrow(new_dat[which(new_dat$group == 1),])),
#                rep(2, nrow(new_dat[which(new_dat$group == 2),])))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#---------------------------------------

#---------------------------------------
# For up genes, for 2 groups

# Restore it under a different name
my_data <- readRDS(paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))
# > my_data[1:4,1:5]
# Cancer_type   Type        ALB     ATP7A       CASP3
# TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519  1.10979667
# TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239 -0.67663103
# TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562  0.01558295
# TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199  0.79229205
# > nrow(my_data)
# [1] 453
# > ncol(my_data)
# [1] 15
my_data <- my_data[which(my_data$Type == "Tumour"),]
my_data$patient_id <- rownames(my_data)
my_data <- my_data %>% dplyr::select(patient_id, everything())
# > my_data[1:4,1:5]
# patient_id Cancer_type   Type        ALB     ATP7A
# TCGA.S9.A7J2.01 TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519
# TCGA.DU.7302.01 TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239
# TCGA.P5.A5EX.01 TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562
# TCGA.HT.8563.01 TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199
# > nrow(my_data)
# [1] 348

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
load(paste(outDir, "/", cancertype, "_surUp.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "1"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "1"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "1"
# [6,] "TCGA.CS.5390.01" "1"

new_dat <- merge(my_data, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
new_dat <- new_dat %>% dplyr::select(group, everything())
new_dat <- new_dat[order(new_dat$group),]
# > new_dat[1:4,1:5]
# group      patient_id Cancer_type   Type        ALB
# 1     1 TCGA.CS.4938.01         LGG Tumour -0.2547535
# 3     1 TCGA.CS.4942.01         LGG Tumour -1.2115040
# 5     1 TCGA.CS.4944.01         LGG Tumour -2.2555637
# 6     1 TCGA.CS.5390.01         LGG Tumour -0.9262113
# > nrow(new_dat)
# [1] 348
# > ncol(new_dat)
# [1] 17

# format the data
dat <- new_dat
# dat <- my_data
dat <- dat[,-c(1,2,3,4)]
colnames(dat) <- gsub("[.]","-",colnames(dat))
dat <- dat[, which(colnames(dat) %in% upGeneList)]
dat <- t(dat)
colnames(dat) <- paste(new_dat[,2], "_", new_dat[,1], sep = "")
# > dat[1:5,1:4]
# TCGA.CS.4938.01_1 TCGA.CS.4942.01_1 TCGA.CS.4944.01_1 TCGA.CS.5390.01_1
# ATP7A       -0.35641294         0.6141353        -0.9737546       0.764118396
# CASP3        0.02466847         0.2716916        -0.5768917      -0.009506082
# CDK1        -0.57554250         0.6503801        -1.7241047       1.358908879
# CP           1.07210424         0.4782524         1.7397376      -0.020783054
# F5           0.77685844         1.5402846         2.3243897       0.669144920
# > nrow(dat)
# [1] 8
# > ncol(dat)
# [1] 348

# For up genes
r <- dat
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% upGeneList),]
#rr <- rbind(rr, r[which(row.names(r) %in% geneList12),])
nGroup1 <- nrow(patients[which(patients[,2] == 1),])
nGroup2 <- nrow(patients[which(patients[,2] == 2),])
# draw
f <- paste(outDir, "/", cancertype, "_heatmapUpGenes.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, nrow(dat)))
# my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, nGroup1), rep(2,nGroup2))
# my_group2 <- c(rep(1, nrow(new_dat[which(new_dat$group == 1),])),
#                rep(2, nrow(new_dat[which(new_dat$group == 2),])))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

#---------------------------------------

#---------------------------------------
# For down genes, for 2 groups

# Restore it under a different name
my_data <- readRDS(paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))
# > my_data[1:4,1:5]
# Cancer_type   Type        ALB     ATP7A       CASP3
# TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519  1.10979667
# TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239 -0.67663103
# TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562  0.01558295
# TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199  0.79229205
# > nrow(my_data)
# [1] 453
# > ncol(my_data)
# [1] 15
my_data <- my_data[which(my_data$Type == "Tumour"),]
my_data$patient_id <- rownames(my_data)
my_data <- my_data %>% dplyr::select(patient_id, everything())
# > my_data[1:4,1:5]
# patient_id Cancer_type   Type        ALB     ATP7A
# TCGA.S9.A7J2.01 TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519
# TCGA.DU.7302.01 TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239
# TCGA.P5.A5EX.01 TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562
# TCGA.HT.8563.01 TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199
# > nrow(my_data)
# [1] 348

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
load(paste(outDir, "/", cancertype, "_surDown.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "2"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "2"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "2"
# [6,] "TCGA.CS.5390.01" "2"

new_dat <- merge(my_data, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
new_dat <- new_dat %>% dplyr::select(group, everything())
new_dat <- new_dat[order(new_dat$group),]
# > new_dat[1:4,1:5]
# group      patient_id Cancer_type   Type         ALB
# 10     1 TCGA.CS.5396.01         LGG Tumour -1.51012541
# 11     1 TCGA.CS.5397.01         LGG Tumour -1.30781119
# 17     1 TCGA.CS.6669.01         LGG Tumour -0.73239273
# 18     1 TCGA.CS.6670.01         LGG Tumour  0.05036363
# > nrow(new_dat)
# [1] 348
# > ncol(new_dat)
# [1] 17

# format the data
dat <- new_dat
# dat <- my_data
dat <- dat[,-c(1,2,3,4)]
colnames(dat) <- gsub("[.]","-",colnames(dat))
dat <- dat[, which(colnames(dat) %in% downGeneList)]
dat <- t(dat)
colnames(dat) <- paste(new_dat[,2], "_", new_dat[,1], sep = "")
# > dat[1:5,1:4]
# TCGA.CS.5396.01_1 TCGA.CS.5397.01_1 TCGA.CS.6669.01_1
# ALB            -1.51012541        -1.3078112        -0.7323927
# CYP1A1          0.82362199         0.3998401         1.4873848
# MAP1LC3A       -1.55906945         0.6137549         1.3166676
# MT-CO1         -1.16969309        -0.5808168         1.3172970
# SNCA            0.01020726         0.7913414         1.7795062
# TCGA.CS.6670.01_1
# ALB             0.05036363
# CYP1A1          0.62906979
# MAP1LC3A        0.09987831
# MT-CO1          0.75338201
# SNCA            0.26303082
# > nrow(dat)
# [1] 5
# > ncol(dat)
# [1] 348

# For down genes
r <- dat
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% downGeneList),]
#rr <- rbind(rr, r[which(row.names(r) %in% geneList12),])
nGroup1 <- nrow(patients[which(patients[,2] == 1),])
nGroup2 <- nrow(patients[which(patients[,2] == 2),])
# draw
f <- paste(outDir, "/", cancertype, "_heatmapDownGenes.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, nrow(dat)))
# my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, nGroup1), rep(2,nGroup2))
# my_group2 <- c(rep(1, nrow(new_dat[which(new_dat$group == 1),])),
#                rep(2, nrow(new_dat[which(new_dat$group == 2),])))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#---------------------------------------

#---------------------------------------
# For up & down genes, for 2 groups

# Restore it under a different name
my_data <- readRDS(paste(outDir, "/", cancertype, "_dataForHeatmap.rds", sep = ""))
# > my_data[1:4,1:5]
# Cancer_type   Type        ALB     ATP7A       CASP3
# TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519  1.10979667
# TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239 -0.67663103
# TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562  0.01558295
# TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199  0.79229205
# > nrow(my_data)
# [1] 453
# > ncol(my_data)
# [1] 15
my_data <- my_data[which(my_data$Type == "Tumour"),]
my_data$patient_id <- rownames(my_data)
my_data <- my_data %>% dplyr::select(patient_id, everything())
# > my_data[1:4,1:5]
# patient_id Cancer_type   Type        ALB     ATP7A
# TCGA.S9.A7J2.01 TCGA.S9.A7J2.01         LGG Tumour  0.2196381 0.2057519
# TCGA.DU.7302.01 TCGA.DU.7302.01         LGG Tumour -0.4343060 0.8561239
# TCGA.P5.A5EX.01 TCGA.P5.A5EX.01         LGG Tumour  0.4626852 0.7158562
# TCGA.HT.8563.01 TCGA.HT.8563.01         LGG Tumour -0.8185620 2.7716199
# > nrow(my_data)
# [1] 348

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
load(paste(outDir, "/", cancertype, "_sur.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "1"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "1"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "1"
# [6,] "TCGA.CS.5390.01" "2"

new_dat <- merge(my_data, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
new_dat <- new_dat %>% dplyr::select(group, everything())
new_dat <- new_dat[order(new_dat$group),]
# > new_dat[1:4,1:5]
# group      patient_id Cancer_type   Type        ALB
# 1      1 TCGA.CS.4938.01         LGG Tumour -0.2547535
# 3      1 TCGA.CS.4942.01         LGG Tumour -1.2115040
# 5      1 TCGA.CS.4944.01         LGG Tumour -2.2555637
# 10     1 TCGA.CS.5396.01         LGG Tumour -1.5101254
# > nrow(new_dat)
# [1] 348
# > ncol(new_dat)
# [1] 17

# format the data
dat <- new_dat
# dat <- my_data
dat <- dat[,-c(1,2,3,4)]
colnames(dat) <- gsub("[.]","-",colnames(dat))
# dat <- dat[, which(colnames(dat) %in% downGeneList)]
dat <- t(dat)
colnames(dat) <- paste(new_dat[,2], "_", new_dat[,1], sep = "")
# > dat[1:5,1:4]
# TCGA.CS.4938.01_1 TCGA.CS.4942.01_1 TCGA.CS.4944.01_1 TCGA.CS.5396.01_1
# ALB         -0.25475346        -1.2115040        -2.2555637       -1.51012541
# ATP7A       -0.35641294         0.6141353        -0.9737546        1.56190764
# CASP3        0.02466847         0.2716916        -0.5768917       -0.21182110
# CDK1        -0.57554250         0.6503801        -1.7241047        0.69373930
# CP           1.07210424         0.4782524         1.7397376        0.02692021
# > nrow(dat)
# [1] 13
# > ncol(dat)
# [1] 348

# For up & down genes
r <- dat
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% upGeneList),]
rr <- rbind(rr, r[which(row.names(r) %in% downGeneList),])
nGroup1 <- nrow(patients[which(patients[,2] == 1),])
nGroup2 <- nrow(patients[which(patients[,2] == 2),])
# draw
f <- paste(outDir, "/", cancertype, "_heatmapUpDownGenes.pdf", sep = "")
pdf(file = f, width = 10, height =  10)
data <- as.matrix(rr)
my_group <- c(rep(1, length(upGeneList)), rep(2, length(downGeneList)))
# my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, nGroup1), rep(2,nGroup2))
# my_group2 <- c(rep(1, nrow(new_dat[which(new_dat$group == 1),])),
#                rep(2, nrow(new_dat[which(new_dat$group == 2),])))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#---------------------------------------

#================================================================

#================================================================
# (25) Umap for tumour classification, using 18 genes and 12 genes for each cancer type, and pan-cancer
# Run in server
#================================================================

# Function to plot umap
plot.umap <- function(x, labels,
         main="A UMAP visualization",
         colors=c("#ff7f00", "#e377c2", "#17becf"),
         pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
         cex.main=1, cex.legend=0.85) {

  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }

  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomleft"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }

  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

# load the data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]
geneList18_12 <- c(geneList18, geneList12)

#---------------------------------------
# For 18 genes, for 2 groups, LGG

cancertype <- "LGG"

outDir <- paste(rootDir, "/Data/Output/surCancerType1812", sep = "")
load(paste(outDir, "/", cancertype, "_sur18.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "2"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "2"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "2"
# [6,] "TCGA.CS.5390.01" "2"

# gene data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

surdata <- surdata[which(surdata$cancertype == cancertype),]
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.CS.4938.01        LGG  0    3574 3.3219
# TCGA.CS.4941.01        LGG  1     234 4.7549
# TCGA.CS.4942.01        LGG  1    1335 3.5850
# TCGA.CS.4943.01        LGG  1    1106 2.0000
# TCGA.CS.4944.01        LGG  0    1828 3.3219
# > nrow(surdata)
# [1] 348
# > ncol(surdata)
# [1] 60

# draw
f <- paste(outDir, "/surCancerType1812/", cancertype, "_umap18genes.pdf", sep = "")
pdf(file = f, width = 5, height =  5)
x <- surdata
x$patient_id <- rownames(x)
x <- merge(x, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
# > x[1:4,1:5]
# patient_id cancertype OS OS.time  AANAT
# 1 TCGA.CS.4938.01        LGG  0    3574 3.3219
# 2 TCGA.CS.4941.01        LGG  1     234 4.7549
# 3 TCGA.CS.4942.01        LGG  1    1335 3.5850
# 4 TCGA.CS.4943.01        LGG  1    1106 2.0000
# > nrow(x)
# [1] 348
# > ncol(x)
# [1] 62
x <- x %>% dplyr::select(group, everything())
dat <- x[,which(colnames(x) %in% geneList18)]
# > dat[1:4,1:5]
# ADAM10   AP1S1     APP    ARF1 ATP6AP1
# 1 11.4615 10.9679 15.9837 13.8728 12.7946
# 2 12.1633 11.8146 16.4526 14.6035 13.5169
# 3 12.2938 12.0532 16.4958 14.1779 13.4029
# 4 12.4757 11.6974 15.9645 14.2748 13.3393
# > nrow(dat)
# [1] 348
# > ncol(dat)
# [1] 18
labels <- x$group
# labels <- paste("Subtype ", labels, sep = "")
head(labels)
length(labels)
# > head(labels)
# [1] "2" "2" "2" "2" "2" "2"
# > length(labels)
# [1] 348
data <- as.matrix(dat)
dat.umap <- umap(data)
plot.umap(dat.umap, labels, main="A UMAP visualisation of subtypes clustered using the 18 genes",
          colors=c("#046C9A","#D69C4E"))
dev.off()

#---------------------------------------

#---------------------------------------
# For 12 genes, for 2 groups, LGG

cancertype <- "LGG"

outDir <- paste(rootDir, "/Data/Output/surCancerType1812", sep = "")
load(paste(outDir, "/", cancertype, "_sur12.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group

# gene data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata

surdata <- surdata[which(surdata$cancertype == cancertype),]

# draw
f <- paste(outDir, "/surCancerType1812/", cancertype, "_umap12genes.pdf", sep = "")
pdf(file = f, width = 5, height =  5)
x <- surdata
x$patient_id <- rownames(x)
x <- merge(x, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
x <- x %>% dplyr::select(group, everything())
dat <- x[,which(colnames(x) %in% geneList12)]
dat[1:4,1:5]
nrow(dat)
ncol(dat)
# > dat[1:4,1:5]
# AOC3    AQP1 CYP1A1 DBH     JUN
# 1 5.5850 14.5588 3.3219   0 12.8684
# 2 5.4263 18.1242 2.0000   3 13.9323
# 3 4.8580 15.1665 6.4757   3 12.1267
# 4 6.4746 12.8332 1.0000   2 13.3166
# > nrow(dat)
# [1] 348
# > ncol(dat)
# [1] 12
labels <- x$group
# labels <- paste("Subtype ", labels, sep = "")
head(labels)
length(labels)
# > head(labels)
# [1] "2" "1" "2" "2" "2" "2"
# > length(labels)
# [1] 348
data <- as.matrix(dat)
dat.umap <- umap(data)
plot.umap(dat.umap, labels, main="A UMAP visualisation of subtypes clustered using the 12 genes",
          colors=c("#046C9A","#D69C4E"))
dev.off()

#---------------------------------------

#---------------------------------------
# For 18 genes, for 2 groups, pan-cancer

outDir <- paste(rootDir, "/Data/Output", sep = "")
load(paste(outDir, "/sur18.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.OR.A5J1.01" "1"
# [2,] "TCGA.OR.A5J2.01" "1"
# [3,] "TCGA.OR.A5J3.01" "1"
# [4,] "TCGA.OR.A5J5.01" "1"
# [5,] "TCGA.OR.A5J6.01" "1"
# [6,] "TCGA.OR.A5J7.01" "1"
# > nrow(patients)
# [1] 6727

# gene data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

# draw
f <- paste(outDir, "/umap18genes.pdf", sep = "")
pdf(file = f, width = 5, height =  5)
x <- surdata
x$patient_id <- rownames(x)
x <- merge(x, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
# > x[1:4,1:5]
# patient_id cancertype OS OS.time  AANAT
# 1 TCGA.02.0047.01        GBM  1     448 4.0000
# 2 TCGA.02.0055.01        GBM  1      76 3.7004
# 3 TCGA.02.2483.01        GBM  0     466 3.4594
# 4 TCGA.02.2485.01        GBM  0     470 2.8074
# > nrow(x)
# [1] 6727
# > ncol(x)
# [1] 62
x <- x %>% dplyr::select(group, everything())
dat <- x[,which(colnames(x) %in% geneList18)]
# > dat[1:4,1:5]
# ADAM10   AP1S1     APP    ARF1 ATP6AP1
# 1 11.8017 10.9440 15.4466 13.5590 12.7785
# 2 12.0382 11.3219 15.5397 14.6902 13.0848
# 3 11.9363 12.1039 16.0381 14.5960 13.5075
# 4 11.7627 11.4041 16.0160 13.7551 13.0755
# > nrow(dat)
# [1] 6727
# > ncol(dat)
# [1] 18
labels <- x$group
# labels <- paste("Subtype ", labels, sep = "")
head(labels)
length(labels)
# > head(labels)
# [1] "1" "1" "1" "1" "1" "1"
# > length(labels)
# [1] 6727
data <- as.matrix(dat)
dat.umap <- umap(data)
plot.umap(dat.umap, labels, main="A UMAP visualisation of subtypes clustered using the 18 genes",
          colors=c("#046C9A","#D69C4E"))
dev.off()

#---------------------------------------

#---------------------------------------
# For 12 genes, for 2 groups, pan-cancer

outDir <- paste(rootDir, "/Data/Output", sep = "")
load(paste(outDir, "/sur12.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.OR.A5J1.01" "1"
# [2,] "TCGA.OR.A5J2.01" "1"
# [3,] "TCGA.OR.A5J3.01" "1"
# [4,] "TCGA.OR.A5J5.01" "1"
# [5,] "TCGA.OR.A5J6.01" "1"
# [6,] "TCGA.OR.A5J7.01" "1"
# > nrow(patients)
# [1] 6727

# gene data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

# draw
f <- paste(outDir, "/umap12genes.pdf", sep = "")
pdf(file = f, width = 5, height =  5)
x <- surdata
x$patient_id <- rownames(x)
x <- merge(x, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
# > x[1:4,1:5]
# patient_id cancertype OS OS.time  AANAT
# 1 TCGA.02.0047.01        GBM  1     448 4.0000
# 2 TCGA.02.0055.01        GBM  1      76 3.7004
# 3 TCGA.02.2483.01        GBM  0     466 3.4594
# 4 TCGA.02.2485.01        GBM  0     470 2.8074
# > nrow(x)
# [1] 6727
# > ncol(x)
# [1] 62
x <- x %>% dplyr::select(group, everything())
dat <- x[,which(colnames(x) %in% geneList12)]
# > dat[1:4,1:5]
# AOC3    AQP1 CYP1A1    DBH     JUN
# 1 7.4429 16.1987  1.000 3.5850 12.5985
# 2 7.7071 11.2096  5.000 3.4594 14.0844
# 3 7.1799 12.7667  3.000 3.8074 14.2015
# 4 6.1293 16.6494  1.585 3.8074 14.9370
# > nrow(dat)
# [1] 6727
# > ncol(dat)
# [1] 12
labels <- x$group
# labels <- paste("Subtype ", labels, sep = "")
head(labels)
length(labels)
# > head(labels)
# [1] "1" "1" "1" "1" "1" "1"
# > length(labels)
# [1] 6727
data <- as.matrix(dat)
dat.umap <- umap(data)
plot.umap(dat.umap, labels, main="A UMAP visualisation of subtypes clustered using the 12 genes",
          colors=c("#046C9A","#D69C4E"))
dev.off()

#---------------------------------------

#================================================================

#================================================================
# (26) Umap for tumour classification, using copper genes for each cancer type
# Run in server
#================================================================

# Function to plot umap
plot.umap <- function(x, labels,
                      main="A UMAP visualization",
                      colors=c("#ff7f00", "#e377c2", "#17becf"),
                      pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                      cex.main=1, cex.legend=0.85) {
  
  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }
  
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomleft"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

cancertype <- "LGG"

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
x <- read.csv(fileName)
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")

# Critical copper genes
criticalCopperGenes <- y
upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'UP'),]$gene
downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'DOWN'),]$gene
geneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] %in% c("UP", "DOWN")),]$gene

#---------------------------------------
# For up genes, for 2 groups

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
load(paste(outDir, "/", cancertype, "_surUp.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "1"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "1"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "1"
# [6,] "TCGA.CS.5390.01" "1"

# gene data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

surdata <- surdata[which(surdata$cancertype == cancertype),]
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.CS.4938.01        LGG  0    3574 3.3219
# TCGA.CS.4941.01        LGG  1     234 4.7549
# TCGA.CS.4942.01        LGG  1    1335 3.5850
# TCGA.CS.4943.01        LGG  1    1106 2.0000
# TCGA.CS.4944.01        LGG  0    1828 3.3219
# > nrow(surdata)
# [1] 348
# > ncol(surdata)
# [1] 60

# draw
f <- paste(outDir, "/surCancerTypeUseCopperGenes/", cancertype, "_umapUpgenes.pdf", sep = "")
pdf(file = f, width = 5, height =  5)
x <- surdata
x$patient_id <- rownames(x)
x <- merge(x, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)
# > x[1:4,1:5]
# patient_id cancertype OS OS.time  AANAT
# 1 TCGA.CS.4938.01        LGG  0    3574 3.3219
# 2 TCGA.CS.4941.01        LGG  1     234 4.7549
# 3 TCGA.CS.4942.01        LGG  1    1335 3.5850
# 4 TCGA.CS.4943.01        LGG  1    1106 2.0000
# > nrow(x)
# [1] 348
# > ncol(x)
# [1] 62
x <- x %>% dplyr::select(group, everything())
dat <- x[,which(colnames(x) %in% upGeneList)]
dat[1:4,1:5]
nrow(dat)
ncol(dat)
# > dat[1:4,1:5]
# ATP7A   CASP3    CDK1      CP      F5
# 1 8.4838 10.3608  6.3923 11.5159 10.0715
# 2 9.1396 10.8556  9.0084 10.2872  8.8074
# 3 9.3772 10.6027  8.9425 10.1364 12.2677
# 4 9.5018 10.4798 12.3250  8.8437  9.6742
# > nrow(dat)
# [1] 348
# > ncol(dat)
# [1] 8
labels <- x$group
# labels <- paste("Subtype ", labels, sep = "")
head(labels)
length(labels)
# > head(labels)
# [1] "1" "2" "1" "2" "1" "1"
# > length(labels)
# [1] 348
data <- as.matrix(dat)
dat.umap <- umap(data)
plot.umap(dat.umap, labels, main="A UMAP visualisation of clustered subtypes",
          colors=c("#046C9A","#D69C4E"), pad=0.1)
dev.off()

#---------------------------------------

#---------------------------------------
# For down genes, for 2 groups

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
load(paste(outDir, "/", cancertype, "_surDown.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "2"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "2"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "2"
# [6,] "TCGA.CS.5390.01" "2"

# gene data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata

surdata <- surdata[which(surdata$cancertype == cancertype),]

# draw
f <- paste(outDir, "/surCancerTypeUseCopperGenes/", cancertype, "_umapDowngenes.pdf", sep = "")
pdf(file = f, width = 5, height =  5)
x <- surdata
x$patient_id <- rownames(x)
x <- merge(x, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)

x <- x %>% dplyr::select(group, everything())
dat <- x[,which(colnames(x) %in% downGeneList)]
dat[1:4,1:5]
nrow(dat)
ncol(dat)
# > dat[1:4,1:5]
# ALB CYP1A1 MAP1LC3A  MT-CO1    SNCA
# 1 5.4263 3.3219  10.6839 20.6123  8.7616
# 2 5.7808 2.0000  11.4057 19.1937 10.2193
# 3 4.0875 6.4757  10.4888 19.3822  9.9159
# 4 3.8074 1.0000   8.5038 19.0112  8.6330
# > nrow(dat)
# [1] 348
# > ncol(dat)
# [1] 5
labels <- x$group
# labels <- paste("Subtype ", labels, sep = "")
head(labels)
length(labels)
# > head(labels)
# [1] "2" "2" "2" "2" "2" "2"
# > length(labels)
# [1] 348
data <- as.matrix(dat)
dat.umap <- umap(data)
plot.umap(dat.umap, labels, main="A UMAP visualisation of clustered subtypes",
          colors=c("#046C9A","#D69C4E"), pad=0.1)
dev.off()

#---------------------------------------

#---------------------------------------
# For up & down genes, for 2 groups

outDir <- paste(rootDir, "/Data/Output/surCancerTypeUseCopperGenes", sep = "")
load(paste(outDir, "/", cancertype, "_sur.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.CS.4938.01" "1"
# [2,] "TCGA.CS.4941.01" "2"
# [3,] "TCGA.CS.4942.01" "1"
# [4,] "TCGA.CS.4943.01" "2"
# [5,] "TCGA.CS.4944.01" "1"
# [6,] "TCGA.CS.5390.01" "2"

# gene data
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata

surdata <- surdata[which(surdata$cancertype == cancertype),]

# draw
f <- paste(outDir, "/surCancerTypeUseCopperGenes/", cancertype, "_umapUpDowngenes.pdf", sep = "")
pdf(file = f, width = 5, height =  5)
x <- surdata
x$patient_id <- rownames(x)
x <- merge(x, patients, by.x = "patient_id", by.y = "patient_id", all.x = TRUE)

x <- x %>% dplyr::select(group, everything())
dat <- x[,which(colnames(x) %in% geneList)]
dat[1:4,1:5]
nrow(dat)
ncol(dat)
# > dat[1:4,1:5]
# ALB  ATP7A   CASP3    CDK1      CP
# 1 5.4263 8.4838 10.3608  6.3923 11.5159
# 2 5.7808 9.1396 10.8556  9.0084 10.2872
# 3 4.0875 9.3772 10.6027  8.9425 10.1364
# 4 3.8074 9.5018 10.4798 12.3250  8.8437
# > nrow(dat)
# [1] 348
# > ncol(dat)
# [1] 13
labels <- x$group
# labels <- paste("Subtype ", labels, sep = "")
head(labels)
length(labels)
# > head(labels)
# [1] "1" "2" "1" "2" "1" "2"
# > length(labels)
# [1] 348
data <- as.matrix(dat)
dat.umap <- umap(data)
plot.umap(dat.umap, labels, main="A UMAP visualisation of clustered subtypes",
          colors=c("#046C9A","#D69C4E"), pad=0.1)
dev.off()

#---------------------------------------

#================================================================

#================================================================
# (27) Compare gene expression, copper genes, 18 & 12 genes, 23 cancer types, using count and cpm
# heatmap for clustered groups in pan-cancer
# identify subset of genes for tumour classification & survival analysis
# run on server
#================================================================

# Link the files
# Katana
# cancertype="ACC"
# ln -s /srv/scratch/z3538133/001pancancer/pan/data/$cancertype/gtex.$cancertype.raw.counts.RData /srv/scratch/z3538133/002NetworkAnalysis/Data/$cancertype/gtex.$cancertype.raw.counts.RData

# # Get the necessary data
# getHeatmapData=function(geneList) {
#   subtypes <- PanCancerAtlas_subtypes()
#   subtypes <- unique(subtypes$cancer.type)
#   subtypes <- c(subtypes, "PAAD")
#   subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
#   subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
#   n <- length(subtypes)
#   Cancer_type <- c()
#   Type <- c()
#   Expression <- c()
#   for (i in 1:n) {
#     print(i)
#     cancertype <- subtypes[i]
#     outDir <- paste(rootDir, "/Data/", cancertype, sep = "")
#     load(paste(outDir, "/gtex.", cancertype, ".raw.counts.RData", sep = "")) # return gtex.PAN.raw.counts, including both tumor and normal
#     # # transform
#     # dat <- gtex.PAN.raw.counts
#     # rownames(dat) <- dat$gene
#     # dat <- dat[,-1]
#     # dat <- dat + 1
#     # dat <- log(dat, base = 2)
#     # iGenes <- which(rownames(dat) %in% geneList)
#     # r <- str_detect(colnames(dat), "TCGA")
#     # nTumor <- length(which(r == TRUE))
#     # nNormal <- ncol(dat) - nTumor
#     
#     # compute logCPM
#     dat <- gtex.PAN.raw.counts
#     rownames(dat) <- dat$gene
#     dat <- dat[,-1]
#     # dat <- dat + 1
#     logCPM <- cpm(dat, prior.count=2, log=TRUE)
#     logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
#     dat <- logCPM
#     iGenes <- which(rownames(dat) %in% geneList)
#     r <- str_detect(colnames(dat), "TCGA")
#     nTumor <- length(which(r == TRUE))
#     nNormal <- ncol(dat) - nTumor
#     
#     # # Add * if being critical copper gene
#     # outDir <- paste(rootDir, "/Data/Output", sep = "")
#     # fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
#     # selectedCriticalCopperGenes <- read.csv(fileName)
#     # if(selectedCriticalCopperGenes[which(selectedCriticalCopperGenes$gene == gene),cancertype] == 1){
#     #   cancertype <- paste("*", cancertype, sep = "")
#     # }
#     
#     Cancer_type <- c(Cancer_type, rep(cancertype, nTumor + nNormal))
#     Type <- c(Type, rep(c("Tumour", "Normal"), times = c(nTumor, nNormal)))
#     x <- dat[iGenes,]
#     x <- t(x)
#     # order by gene names
#     x <- x[, order(colnames(x), decreasing = FALSE)]
#     Expression <- rbind(Expression, x)
#     
#   }
#   data=data.frame(Cancer_type, Type ,  Expression)
#   
#   return(data)
# }

outDir <- paste(rootDir, "/Data/Output", sep = "")

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# nrow(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# [1] 35

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]

# 30 genes
# list30 <- c(geneList18, geneList12)
# cat(list30, sep = ",")
# CDK1,AP1S1,CASP3,TMPRSS6,GSK3B,APP,COX17,XIAP,ARF1,GPC1,SORD,ATP7A,SP1,MT-CO1,SLC11A2,ATP6AP1,ADAM10,CP,MAP1LC3A,SNCA,MAPT,JUN,CYP1A1,AOC3,PRNP,DBH,S100A12,AQP1,MT1X,XAF1
# > cat(geneList18, sep = ",")
# CDK1,AP1S1,CASP3,TMPRSS6,GSK3B,APP,COX17,XIAP,ARF1,GPC1,SORD,ATP7A,SP1,MT-CO1,SLC11A2,ATP6AP1,ADAM10,CP
#   > cat(geneList12, sep = ",")
# MAP1LC3A,SNCA,MAPT,JUN,CYP1A1,AOC3,PRNP,DBH,S100A12,AQP1,MT1X,XAF1

# Restore it under a different name
my_data <- readRDS(paste(outDir, "/dataForHeatmap.rds", sep = ""))
# > my_data[1:5, 1:6]
# Cancer_type   Type     ADAM10       AOC3      AP1S1        APP
# TCGA.OR.A5K3.01         ACC Tumour -0.7909516 -2.8371850 -0.2222094 -1.1530789
# TCGA.OR.A5J2.01         ACC Tumour  1.4904894 -0.8957145  1.2599977  0.7926713
# TCGA.OR.A5LN.01         ACC Tumour -0.9779454 -1.9035118 -1.0292136 -0.5072996
# TCGA.OR.A5KY.01         ACC Tumour  0.3102540 -0.3041889  2.3457249  0.6445042
# TCGA.OR.A5LG.01         ACC Tumour  0.7147539  0.5562921  0.6389964 -0.3637465
# > nrow(my_data)
# [1] 11329
# > ncol(my_data)
# [1] 32

# Only get tumour data
my_data <- my_data[which(my_data[,2] == "Tumour"),]
nrow(my_data)
# nrow(my_data)
# [1] 6732

# format the data
dat <- my_data
dat <- dat[,-c(1,2)] # remove columns of cancer type and type of normal or tumour
dat <- t(dat)

r <- dat

# Divide by up genes and down genes
# geneList18 <- gsub("[-]",".",geneList18)
# geneList12 <- gsub("[-]",".",geneList12)
rownames(r) <- gsub("[.]","-",rownames(r))
rr <- r[which(row.names(r) %in% geneList18),]
rr <- rbind(rr,r[which(row.names(r) %in% geneList12),])
# > rr[1:4,1:5]
# TCGA.OR.A5K3.01 TCGA.OR.A5J2.01 TCGA.OR.A5LN.01 TCGA.OR.A5KY.01
# ADAM10      -0.7909516       1.4904894      -0.9779454       0.3102540
# AP1S1       -0.2222094       1.2599977      -1.0292136       2.3457249
# APP         -1.1530789       0.7926713      -0.5072996       0.6445042
# ARF1        -0.5307502       1.1526288      -0.1879886       0.4796551
# TCGA.OR.A5LG.01
# ADAM10       0.7147539
# AP1S1        0.6389964
# APP         -0.3637465
# ARF1         0.5872978
ncol(rr)
# > ncol(rr)
# [1] 6732
write.csv(rr, paste(outDir, "/rr_heatmap_pan_groups.csv", sep = ""), row.names=TRUE)

#---------------------------------------
# For 18 genes, for 2 groups, pan-cancer

outDir <- paste(rootDir, "/Data/Output", sep = "")
load(paste(outDir, "/sur18.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.OR.A5J1.01" "1"
# [2,] "TCGA.OR.A5J2.01" "1"
# [3,] "TCGA.OR.A5J3.01" "1"
# [4,] "TCGA.OR.A5J5.01" "1"
# [5,] "TCGA.OR.A5J6.01" "1"
# [6,] "TCGA.OR.A5J7.01" "1"
# > nrow(patients)
# [1] 6727

# # gene data
# outDir <- paste(rootDir, "/Data/Output", sep = "")
# fileName<-paste(outDir, "/surdata.Rdata",sep="")
# load(fileName) # return surdata
# surdata[1:5,1:4]
# nrow(surdata)
# ncol(surdata)
# # > surdata[1:5,1:4]
# # cancertype OS OS.time  AANAT
# # TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# # TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# # TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# # TCGA.OR.A5J5.01        ACC  1     365 2.8074
# # TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# # > nrow(surdata)
# # [1] 6727
# # > ncol(surdata)
# # [1] 60

rr <- rr[which(row.names(rr) %in% geneList18),]
r1 <- rr[,which(colnames(rr) %in% patients[which(patients[,2] == 1),1])]
r2 <- cbind(r1, rr[,which(colnames(rr) %in% patients[which(patients[,2] == 2),1])])
# > ncol(r1)
# [1] 6537
# > ncol(r2) - ncol(r1)
# [1] 190

# draw the chart
f <- paste(outDir, "/heatmap18_pan.png", sep = "")
png(file = f, width = 600, height =  600)
data <- as.matrix(r2)
my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, 6537), rep(2,190))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

x <- surdata
x <- x[,which(colnames(x) %in% geneList18)]
# change from log count to count
x_count <- 2^x
# > x_count[1:4,1:5]
# ADAM10     AP1S1      APP     ARF1  ATP6AP1
# TCGA.OR.A5J1.01 3600.0469  1621.036 51686.50 24648.29 11853.29
# TCGA.OR.A5J2.01 6566.9147  6434.891 94668.71 48466.49 16760.76
# TCGA.OR.A5J3.01 1924.9430 11879.611 75783.70 18333.57 16463.69
# TCGA.OR.A5J5.01  946.9865  2803.085 23731.24 20122.19  8097.72
# > nrow(x_count)
# [1] 6727
# > ncol(x_count)
# [1] 18

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
dat <- logCPM
# > dat[1:4,1:5]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01 TCGA.OR.A5J6.01
# ADAM10       -2.175977     -0.90739836     -1.98254127       -4.269179      -1.5033557
# AP1S1        -1.896608      0.09502618      1.07297573       -1.814749      -0.3233564
# APP          -1.109462     -0.12276656     -0.05020926       -2.235250      -1.0692236
# ARF1         -1.853613     -0.35860223     -1.21688994       -2.717414      -1.4324985
# > nrow(dat)
# [1] 18
# > ncol(dat)
# [1] 6727

rr <- dat
r1 <- rr[,which(colnames(rr) %in% patients[which(patients[,2] == 1),1])]
r2 <- cbind(r1, rr[,which(colnames(rr) %in% patients[which(patients[,2] == 2),1])])
# > ncol(r1)
# [1] 6537
# > ncol(r2) - ncol(r1)
# [1] 190

# draw the chart
f <- paste(outDir, "/heatmap18_pan.png", sep = "")
png(file = f, width = 600, height =  600)
data <- as.matrix(r2)
my_group <- c(rep(1, 18))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, 6537), rep(2,190))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#-----------------------------------------

#---------------------------------------
# For 12 genes, for 2 groups, pan-cancer

outDir <- paste(rootDir, "/Data/Output", sep = "")

rr <- read.csv(paste(outDir, "/rr_heatmap_pan_groups.csv", sep = ""), row.names = 1)

load(paste(outDir, "/sur12.Rdata", sep = "")) # return result1

patients <- matrix(NA, nrow = nrow(result1$distanceMatrix), ncol = 2)
colnames(patients) <- c("patient_id", "group")
patients[,1] <- as.character(colnames(result1$distanceMatrix))
patients[,2] <- result1$group
# > head(patients)
# patient_id        group
# [1,] "TCGA.OR.A5J1.01" "1"
# [2,] "TCGA.OR.A5J2.01" "1"
# [3,] "TCGA.OR.A5J3.01" "1"
# [4,] "TCGA.OR.A5J5.01" "1"
# [5,] "TCGA.OR.A5J6.01" "1"
# [6,] "TCGA.OR.A5J7.01" "1"
# > nrow(patients)
# [1] 6727

rr <- rr[which(row.names(rr) %in% geneList12),]
r1 <- rr[,which(colnames(rr) %in% patients[which(patients[,2] == 1),1])]
r2 <- cbind(r1, rr[,which(colnames(rr) %in% patients[which(patients[,2] == 2),1])])
# > ncol(r1)
# [1] 6553
# > ncol(r2) - ncol(r1)
# [1] 174

# draw the chart
f <- paste(outDir, "/heatmap12_pan.png", sep = "")
png(file = f, width = 600, height =  600)
data <- as.matrix(r2)
my_group <- c(rep(1, 12))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, 6553), rep(2,174))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
#ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

x <- surdata
x <- x[,which(colnames(x) %in% geneList12)]
# change from log count to count
x_count <- 2^x
# > x_count[1:4,1:5]
# AOC3      AQP1    CYP1A1 DBH      JUN
# TCGA.OR.A5J1.01  260.9989 1575.3973  3.000078   1 5587.937
# TCGA.OR.A5J2.01 1273.9562 2840.6394 17.999688   8 4454.970
# TCGA.OR.A5J3.01  625.8617  440.9893  8.999844   4 4414.089
# TCGA.OR.A5J5.01  528.6247  576.9890 16.000000   1 3019.094
# > nrow(x_count)
# [1] 6727
# > ncol(x_count)
# [1] 12

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
dat <- logCPM
# > dat[1:4,1:5]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01
# AOC3         -0.665191      0.82017115       0.1657128       0.2340588
# AQP1         -1.000626      0.08404926      -1.8267162      -1.3414037
# CYP1A1       -0.374508      0.89360805       0.3588996       0.8334210
# DBH          -1.199724      0.03423468      -0.4766345      -0.9840322
# TCGA.OR.A5J6.01
# AOC3         1.4112064
# AQP1         0.3413380
# CYP1A1       0.2164131
# DBH         -0.5964999
# > nrow(dat)
# [1] 12
# > ncol(dat)
# [1] 6727

rr <- dat
r1 <- rr[,which(colnames(rr) %in% patients[which(patients[,2] == 1),1])]
r2 <- cbind(r1, rr[,which(colnames(rr) %in% patients[which(patients[,2] == 2),1])])
# > ncol(r1)
# [1] 6553
# > ncol(r2) - ncol(r1)
# [1] 174

# draw the chart
f <- paste(outDir, "/heatmap12_pan.png", sep = "")
png(file = f, width = 600, height =  600)
data <- as.matrix(r2)
my_group <- c(rep(1, 12))
rowSide <- brewer.pal(3, "Set1")[my_group]
my_group2 <- c(rep(1, 6553), rep(2,174))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#-----------------------------------------

#---------------------------------------
#------------------------------------
# For 2 groups only
#------------------------------------
# For 30 genes, identify 2 groups up & down, pan-cancer
# Top5, including both up and down

outDir <- paste(rootDir, "/Data/Output", sep = "")

# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the cases of 30 genes, including 18 genes and 12 genes, for both up and down genes
exp_data <- surdata[,-c(1,2,3)]
exp_data <- exp_data[,colnames(exp_data) %in% c(geneList18, geneList12)]
exp_data[1:4,1:5]
nrow(exp_data)
ncol(exp_data)
# ADAM10    AOC3   AP1S1     APP    AQP1
# TCGA.OR.A5J1.01 11.8138  8.0279 10.6627 15.6575 10.6215
# TCGA.OR.A5J2.01 12.6810 10.3151 12.6517 16.5306 11.4720
# TCGA.OR.A5J3.01 10.9106  9.2897 13.5362 16.2096  8.7846
# TCGA.OR.A5J5.01  9.8872  9.0461 11.4528 14.5345  9.1724
# [1] 6727
# [1] 30

x <- exp_data

# change from log count to count
x_count <- 2^x
# > x_count[1:4,1:5]
# ADAM10      AOC3     AP1S1      APP      AQP1
# TCGA.OR.A5J1.01 3600.0469  260.9989  1621.036 51686.50 1575.3973
# TCGA.OR.A5J2.01 6566.9147 1273.9562  6434.891 94668.71 2840.6394
# TCGA.OR.A5J3.01 1924.9430  625.8617 11879.611 75783.70  440.9893
# TCGA.OR.A5J5.01  946.9865  528.6247  2803.085 23731.24  576.9890
# > nrow(x_count)
# [1] 6727
# > ncol(x_count)
# [1] 30

# compute CPM
dat <- t(x_count)
CPM <- cpm(dat, log=FALSE)

exp_data <- CPM
exp_data <- t(exp_data)
mean30 <- colMeans(exp_data)
# > mean30
# ADAM10        AOC3       AP1S1         APP        AQP1        ARF1
# 7252.3709   2600.3428   3196.4260  59134.3081  16067.4004  30894.1273
# ATP6AP1       ATP7A       CASP3        CDK1       COX17          CP
# 13193.3335   1309.5462   2062.9998   2080.0265   1825.1825  10430.6465
# CYP1A1         DBH        GPC1       GSK3B         JUN    MAP1LC3A
# 266.7921   4534.4520   6674.0926   4320.8433  14387.6737   1209.0188
# MAPT      MT-CO1        MT1X        PRNP     S100A12     SLC11A2
# 3583.1980 775925.9396   3771.6661   8363.3308    103.0065   4360.4695
# SNCA        SORD         SP1     TMPRSS6        XAF1        XIAP
# 860.8214   7855.0815   7069.0119    702.0836   2211.8583   3753.9499

result_list <- list()
list_name <- c(1:30)
for (i in 1:30) {
  gene <- colnames(exp_data)[i]
  high <- rownames(exp_data)[which(exp_data[,i] >= mean30[i])]
  low <- rownames(exp_data)[which(exp_data[,i] < mean30[i])]
  result_list[[list_name[i]]] <- list(gene = gene, high = high, low = low)
}

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y <- y[which(y$frequency > 1),]
y_plus <- y
y_plus$UpDown <- "Varied"
y_plus$UpDown <- ifelse(y_plus$numUp > y_plus$numDown, "Up", y_plus$UpDown)
y_plus$UpDown <- ifelse(y_plus$numUp < y_plus$numDown, "Down", y_plus$UpDown)
y_plus$percent <- ifelse(y_plus$UpDown == "Up", y_plus$numUp/y_plus$frequency, y_plus$numDown/y_plus$frequency)
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
y_plus <- y_plus[which(y_plus$gene %in% uni_gene$x),]

# gene_set <- c("ATP7A", "CP", "APP", "MAPT", "DBH", "TMPRSS6", "AOC3", "ADAM10", "SP1", "XIAP")
# gene_set <- c("CDK1", "AP1S1", "CASP3")
# gene_set <- y_plus[which(y_plus$UpDown == "Up" & y_plus$percent >= 1),]$gene
gene_set <- y_plus[which(y_plus$frequency > 7),]$gene # top 5
# > gene_set
# [1] "CDK1"     "AP1S1"    "CASP3"    "MAP1LC3A" "SNCA"      
# > result_list[[3]]$gene
# [1] "AP1S1"
# "CDK1"     "AP1S1"    "CASP3" high & "MAP1LC3A" "SNCA" low
res_high <- result_list[[3]]$high # index 3 for AP1S1
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    if(result_list[[i]]$gene %in% c("CDK1",     "AP1S1",    "CASP3")) {
      res_high <- intersect(res_high, result_list[[i]]$high)  
    } else {
      res_high <- intersect(res_high, result_list[[i]]$low)
    }
    print(length(res_high))
  }
}
# [1] 2283
# [1] 1620
# [1] 1177
# [1] 493
# [1] 390

# "CDK1"     "AP1S1"    "CASP3" low and "MAP1LC3A" "SNCA" high
# > result_list[[18]]$gene
# [1] "MAP1LC3A"
res_low <- result_list[[18]]$high # index 18 for MAP1LC3A
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    if(result_list[[i]]$gene %in% c( "CDK1",     "AP1S1",    "CASP3")) {
      res_low <- intersect(res_low, result_list[[i]]$low)  
    } else {
      res_low <- intersect(res_low, result_list[[i]]$high)
    }
      
    print(length(res_low))
  }
}
# [1] 868
# [1] 655
# [1] 568
# [1] 568
# [1] 138

# tumour classification & survival analysis
surdata$group <- 3
for (i in 1:nrow(surdata)) {
  if(row.names(surdata)[i] %in% res_high) {
    surdata$group[i] <- 1
  }
  if(row.names(surdata)[i] %in% res_low) {
    surdata$group[i] <- 2
  }
}
nrow(surdata[which(surdata$group == 1),]) # 3 genes high 2 genes low
nrow(surdata[which(surdata$group == 2),]) # 3 genes low 2 genes high
nrow(surdata[which(surdata$group == 3),]) # others
# [1] 390
# [1] 138
# [1] 6199

surdata <- surdata[which(surdata$group %in% c(1,2)),]
clusterNum = length(unique(surdata$group))
time <- surdata$OS.time
status <- surdata$OS
group <- surdata$group
dataset = list(time, status, x = group)  
surv = survfit(Surv(time, status) ~ x, dataset)
mainTitle <- "Survival analysis"
if(clusterNum>1){
  sdf=NULL
  sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
  print(sdf)
  p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
}else{
  cat("There is only one cluster in the group")
  p_value=1
}

# survival chart
# draw the chart
f <- paste(outDir, "/sur_chart_3up2down.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
# myCol <- wes_palette("Zissou1")
myCol <- wes_palette("Darjeeling2")

# Graph 1
mar.default <- c(5,4,4,2) + 0.1 # c(bottom, left, top, right)
par(mar = mar.default + c(0, 1, 0, 0)) 
title=paste(mainTitle, " (", clusterNum, " clusters)", sep="")
plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
     main = title, font.main=2, cex.lab=1.6, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
       lty=1, lwd=3, cex=1.3, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
       seg.len = 0.3)
digit=ceiling(-log10(p_value)+2)
if (p_value < 2e-16) {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value < 2e-16"),col="blue",font=2,cex=1.5)
} else {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1.5)  
}
dev.off()

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

x <- surdata
x <- x[,which(colnames(x) %in% c(geneList18, geneList12))]

# change from log count to count
x_count <- 2^x
# > x_count[1:4,1:5]
# ADAM10      AOC3     AP1S1      APP      AQP1
# TCGA.OR.A5J1.01 3600.0469  260.9989  1621.036 51686.50 1575.3973
# TCGA.OR.A5J2.01 6566.9147 1273.9562  6434.891 94668.71 2840.6394
# TCGA.OR.A5J3.01 1924.9430  625.8617 11879.611 75783.70  440.9893
# TCGA.OR.A5J5.01  946.9865  528.6247  2803.085 23731.24  576.9890
# > nrow(x_count)
# [1] 6727
# > ncol(x_count)
# [1] 30

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
dat <- logCPM
# > dat[1:4,1:5]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01
# ADAM10       -2.146292     -0.85689015    -1.955022627       -4.266049
# AOC3         -1.626252     -0.33319563    -0.592047968       -1.478044
# AP1S1        -1.875117      0.15632713     1.147479061       -1.785834
# APP          -1.084996     -0.07517397    -0.005885496       -2.231244
# TCGA.OR.A5J6.01
# ADAM10      -1.4791372
# AOC3         0.9058438
# AP1S1       -0.2858717
# APP         -1.0562735
# > nrow(dat)
# [1] 30
# > ncol(dat)
# [1] 6727

rr <- dat
rr <- rr[which(row.names(rr) %in% gene_set),]
r1 <- rr[,which(colnames(rr) %in% res_high)]
r2 <- cbind(r1, rr[,which(colnames(rr) %in% res_low)])
# > ncol(r1)
# [1] 390
# > ncol(r2) - ncol(r1)
# [1] 138

# draw the chart
f <- paste(outDir, "/heatmap_3up2down.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
data <- as.matrix(r2)
my_group <- c(rep(1, 5))
rowSide <- brewer.pal(3, "Set1")[my_group]
#my_group2 <- c(rep(1, 489), rep(2,390))
my_group2 <- c(rep(1, 390), rep(2,138))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#------------------------------------

#------------------------------------
# For 3 groups
#------------------------------------
# For 30 genes, identify 2 groups up & down, pan-cancer
# Top5, including both up and down

outDir <- paste(rootDir, "/Data/Output", sep = "")

# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the cases of 30 genes, including 18 genes and 12 genes, for both up and down genes
exp_data <- surdata[,-c(1,2,3)]
exp_data <- exp_data[,colnames(exp_data) %in% c(geneList18, geneList12)]
exp_data[1:4,1:5]
nrow(exp_data)
ncol(exp_data)
# ADAM10    AOC3   AP1S1     APP    AQP1
# TCGA.OR.A5J1.01 11.8138  8.0279 10.6627 15.6575 10.6215
# TCGA.OR.A5J2.01 12.6810 10.3151 12.6517 16.5306 11.4720
# TCGA.OR.A5J3.01 10.9106  9.2897 13.5362 16.2096  8.7846
# TCGA.OR.A5J5.01  9.8872  9.0461 11.4528 14.5345  9.1724
# [1] 6727
# [1] 30

x <- exp_data

# change from log count to count
x_count <- 2^x

# compute CPM
dat <- t(x_count)
CPM <- cpm(dat, log=FALSE)

exp_data <- CPM
exp_data <- t(exp_data)
mean30 <- colMeans(exp_data)

result_list <- list()
list_name <- c(1:30)
for (i in 1:30) {
  gene <- colnames(exp_data)[i]
  high <- rownames(exp_data)[which(exp_data[,i] >= mean30[i])]
  low <- rownames(exp_data)[which(exp_data[,i] < mean30[i])]
  result_list[[list_name[i]]] <- list(gene = gene, high = high, low = low)
}

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y <- y[which(y$frequency > 1),]
y_plus <- y
y_plus$UpDown <- "Varied"
y_plus$UpDown <- ifelse(y_plus$numUp > y_plus$numDown, "Up", y_plus$UpDown)
y_plus$UpDown <- ifelse(y_plus$numUp < y_plus$numDown, "Down", y_plus$UpDown)
y_plus$percent <- ifelse(y_plus$UpDown == "Up", y_plus$numUp/y_plus$frequency, y_plus$numDown/y_plus$frequency)
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
y_plus <- y_plus[which(y_plus$gene %in% uni_gene$x),]

# gene_set <- c("ATP7A", "CP", "APP", "MAPT", "DBH", "TMPRSS6", "AOC3", "ADAM10", "SP1", "XIAP")
# gene_set <- c("CDK1", "AP1S1", "CASP3")
# gene_set <- y_plus[which(y_plus$UpDown == "Up" & y_plus$percent >= 1),]$gene
gene_set <- y_plus[which(y_plus$frequency > 7),]$gene # top 5
# > gene_set
# [1] "CDK1"     "AP1S1"    "CASP3"    "MAP1LC3A" "SNCA"      
# > result_list[[3]]$gene
# [1] "AP1S1"
# "CDK1"     "AP1S1"    "CASP3" high & "MAP1LC3A" "SNCA" low
res_high <- result_list[[3]]$high # index 3 for AP1S1
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    if(result_list[[i]]$gene %in% c("CDK1",     "AP1S1",    "CASP3")) {
      res_high <- intersect(res_high, result_list[[i]]$high)  
    } else {
      res_high <- intersect(res_high, result_list[[i]]$low)
    }
    print(length(res_high))
  }
}
# [1] 2283
# [1] 1620
# [1] 1177
# [1] 493
# [1] 390

# "CDK1"     "AP1S1"    "CASP3" low and "MAP1LC3A" "SNCA" high
# > result_list[[18]]$gene
# [1] "MAP1LC3A"
res_low <- result_list[[18]]$high # index 18 for MAP1LC3A
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    if(result_list[[i]]$gene %in% c( "CDK1",     "AP1S1",    "CASP3")) {
      res_low <- intersect(res_low, result_list[[i]]$low)  
    } else {
      res_low <- intersect(res_low, result_list[[i]]$high)
    }
    
    print(length(res_low))
  }
}
# [1] 868
# [1] 655
# [1] 568
# [1] 568
# [1] 138

# tumour classification & survival analysis
surdata$group <- 3
for (i in 1:nrow(surdata)) {
  if(row.names(surdata)[i] %in% res_high) {
    surdata$group[i] <- 1
  }
  if(row.names(surdata)[i] %in% res_low) {
    surdata$group[i] <- 2
  }
}
nrow(surdata[which(surdata$group == 1),]) # 3 genes high 2 genes low
nrow(surdata[which(surdata$group == 2),]) # 3 genes low 2 genes high
nrow(surdata[which(surdata$group == 3),]) # others
# [1] 390
# [1] 138
# [1] 6199

# surdata <- surdata[which(surdata$group %in% c(1,2)),]
clusterNum = length(unique(surdata$group))
time <- surdata$OS.time
status <- surdata$OS
group <- surdata$group
dataset = list(time, status, x = group)  
surv = survfit(Surv(time, status) ~ x, dataset)
mainTitle <- "Survival analysis"
if(clusterNum>1){
  sdf=NULL
  sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
  print(sdf)
  p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
}else{
  cat("There is only one cluster in the group")
  p_value=1
}

# survival chart
# draw the chart
f <- paste(outDir, "/sur_chart_3up2down_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
# myCol <- wes_palette("Zissou1")
myCol <- wes_palette("Darjeeling2")

# Graph 1
mar.default <- c(5,4,4,2) + 0.1 # c(bottom, left, top, right)
par(mar = mar.default + c(0, 1, 0, 0)) 
title=paste(mainTitle, " (", clusterNum, " clusters)", sep="")
plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
     main = title, font.main=2, cex.lab=1.6, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
       lty=1, lwd=3, cex=1.3, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
       seg.len = 0.3)
digit=ceiling(-log10(p_value)+2)
if (p_value < 2e-16) {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value < 2e-16"),col="blue",font=2,cex=1.5)
} else {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1.5)  
}
dev.off()

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

x <- surdata
x <- x[,which(colnames(x) %in% c(geneList18, geneList12))]

# change from log count to count
x_count <- 2^x

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
dat <- logCPM

rr <- dat
rr <- rr[which(row.names(rr) %in% gene_set),]
r1 <- rr[,which(colnames(rr) %in% res_high)]
ncol(r1)
r2 <- cbind(r1, rr[,which(colnames(rr) %in% res_low)])
ncol(r2) - ncol(r1)
r3 <- cbind(r2, rr[,which(!(colnames(rr) %in% c(res_high,res_low)))])
ncol(r3) - ncol(r2)
r2 <- r3
# [1] 390
# [1] 138
# [1] 6199

# draw the chart
f <- paste(outDir, "/heatmap_3up2down_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
data <- as.matrix(r2)
my_group <- c(rep(1, 5))
rowSide <- brewer.pal(3, "Set1")[my_group]
#my_group2 <- c(rep(1, 489), rep(2,390))
my_group2 <- c(rep(1, 390), rep(2,138), rep(3,6199))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

#----------------------------
# check cancer types
# group 1
g1 <- surdata[which(row.names(surdata) %in% res_high),]
nrow(g1)
table(g1$cancertype)
# nrow(g1)
# table(g1$cancertype)
# [1] 390
# 
# ACC BLCA BRCA COAD ESCA  GBM KIRC KIRP  LGG LIHC LUAD LUSC OVCA PAAD SKCM STAD
# 1   32  162   73    6    3    2    2    2    6   39   28    6    4    2    7
# UCEC  UCS
# 13    2

# group 2
g2 <- surdata[which(row.names(surdata) %in% res_low),]
nrow(g2)
table(g2$cancertype)
# nrow(g2)
# table(g2$cancertype)
# [1] 138
# 
# GBM KICH  LGG LUAD OVCA PAAD PCPG PRAD SKCM THCA UCEC
# 7    2   56    1    3    1   24   23    9    8    4

# group 3
g3 <- surdata[which(!(colnames(rr) %in% c(res_high,res_low))),]
nrow(g3)
table(g3$cancertype)
# [1] 6199
# 
# ACC  AML BLCA BRCA COAD ESCA  GBM KICH KIRC KIRP  LGG LIHC LUAD LUSC OVCA PAAD
# 75  172  112 1049  204  173  155   89  506  187  290  233  234  164  409  148
# PCPG PRAD SKCM STAD THCA UCEC  UCS
# 155  350  320  402  552  165   55
# [1] 6199

#------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------
# 3 up
#------------------------------------
# the cases of 30 genes, including 18 genes and 12 genes, for up genes in 18 genes
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata

exp_data <- surdata[,-c(1,2,3)]
exp_data <- exp_data[,colnames(exp_data) %in% c(geneList18, geneList12)]
exp_data[1:4,1:5]
nrow(exp_data)
ncol(exp_data)
# > exp_data[1:4,1:5]
# ADAM10    AOC3   AP1S1     APP    AQP1
# TCGA.OR.A5J1.01 11.8138  8.0279 10.6627 15.6575 10.6215
# TCGA.OR.A5J2.01 12.6810 10.3151 12.6517 16.5306 11.4720
# TCGA.OR.A5J3.01 10.9106  9.2897 13.5362 16.2096  8.7846
# TCGA.OR.A5J5.01  9.8872  9.0461 11.4528 14.5345  9.1724
# > nrow(exp_data)
# [1] 6727
# > ncol(exp_data)
# [1] 30

x <- exp_data

# change from log count to count
x_count <- 2^x
# > x_count[1:4,1:5]
# ADAM10      AOC3     AP1S1      APP      AQP1
# TCGA.OR.A5J1.01 3600.0469  260.9989  1621.036 51686.50 1575.3973
# TCGA.OR.A5J2.01 6566.9147 1273.9562  6434.891 94668.71 2840.6394
# TCGA.OR.A5J3.01 1924.9430  625.8617 11879.611 75783.70  440.9893
# TCGA.OR.A5J5.01  946.9865  528.6247  2803.085 23731.24  576.9890
# > nrow(x_count)
# [1] 6727
# > ncol(x_count)
# [1] 30

# compute CPM
dat <- t(x_count)
CPM <- cpm(dat, log=FALSE)

exp_data <- CPM
exp_data <- t(exp_data)
mean30 <- colMeans(exp_data)
# > mean30
# ADAM10        AOC3       AP1S1         APP        AQP1        ARF1     ATP6AP1       ATP7A       CASP3        CDK1       COX17          CP 
# 7252.3709   2600.3428   3196.4260  59134.3081  16067.4004  30894.1273  13193.3335   1309.5462   2062.9998   2080.0265   1825.1825  10430.6465 
# CYP1A1         DBH        GPC1       GSK3B         JUN    MAP1LC3A        MAPT      MT-CO1        MT1X        PRNP     S100A12     SLC11A2 
# 266.7921   4534.4520   6674.0926   4320.8433  14387.6737   1209.0188   3583.1980 775925.9396   3771.6661   8363.3308    103.0065   4360.4695 
# SNCA        SORD         SP1     TMPRSS6        XAF1        XIAP 
# 860.8214   7855.0815   7069.0119    702.0836   2211.8583   3753.9499     

result_list <- list()
list_name <- c(1:30)
for (i in 1:30) {
  gene <- colnames(exp_data)[i]
  high <- rownames(exp_data)[which(exp_data[,i] >= mean30[i])]
  low <- rownames(exp_data)[which(exp_data[,i] < mean30[i])]
  result_list[[list_name[i]]] <- list(gene = gene, high = high, low = low)
}

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y <- y[which(y$frequency > 1),]
y_plus <- y
y_plus$UpDown <- "Varied"
y_plus$UpDown <- ifelse(y_plus$numUp > y_plus$numDown, "Up", y_plus$UpDown)
y_plus$UpDown <- ifelse(y_plus$numUp < y_plus$numDown, "Down", y_plus$UpDown)
y_plus$percent <- ifelse(y_plus$UpDown == "Up", y_plus$numUp/y_plus$frequency, y_plus$numDown/y_plus$frequency)
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
y_plus <- y_plus[which(y_plus$gene %in% uni_gene$x),]

# gene_set <- c("ATP7A", "CP", "APP", "MAPT", "DBH", "TMPRSS6", "AOC3", "ADAM10", "SP1", "XIAP")
# gene_set <- c("CDK1", "AP1S1", "CASP3")
# gene_set <- y_plus[which(y_plus$UpDown == "Up" & y_plus$percent >= 1),]$gene
# gene_set <- y_plus[which(y_plus$frequency > 7),]$gene # top 5
gene_set <- y_plus[which(y_plus$frequency > 7 & y_plus$UpDown == "Up"),]$gene
gene_set
result_list[[3]]$gene
# > gene_set
# [1] "CDK1"  "AP1S1" "CASP3"
# > result_list[[3]]$gene
# [1] "AP1S1"
# all genes up
res_high <- result_list[[3]]$high
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    res_high <- intersect(res_high, result_list[[i]]$high)  
    print(length(res_high))
  }
}
# [1] 2283
# [1] 1620
# [1] 1177

# genes low
# > result_list[[3]]$gene
# [1] "AP1S1"
res_low <- result_list[[3]]$low # index 3 for AP1S1
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    res_low <- intersect(res_low, result_list[[i]]$low)  
    print(length(res_low))
  }
}
# [1] 4444
# [1] 3637
# [1] 3296

# tumour classification & survival analysis
surdata$group <- 3
for (i in 1:nrow(surdata)) {
  if(row.names(surdata)[i] %in% res_high) {
    surdata$group[i] <- 1
  }
  if(row.names(surdata)[i] %in% res_low) {
    surdata$group[i] <- 2
  }
}
nrow(surdata[which(surdata$group == 1),]) # 3 genes up
nrow(surdata[which(surdata$group == 2),]) # 3 genes down
# > nrow(surdata[which(surdata$group == 1),]) # 3 genes up
# [1] 1177
# > nrow(surdata[which(surdata$group == 2),]) # 3 genes down
# [1] 3296
surdata <- surdata[which(surdata$group %in% c(1,2)),]
clusterNum = length(unique(surdata$group))
time <- surdata$OS.time
status <- surdata$OS
group <- surdata$group
dataset = list(time, status, x = group)  
surv = survfit(Surv(time, status) ~ x, dataset)
mainTitle <- "Survival analysis"
if(clusterNum>1){
  sdf=NULL
  sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
  print(sdf)
  p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
}else{
  cat("There is only one cluster in the group")
  p_value=1
}

# survival chart
# draw the chart
f <- paste(outDir, "/sur_chart_3up.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
# myCol <- wes_palette("Zissou1")
myCol <- wes_palette("Darjeeling2")

# Graph 1
mar.default <- c(5,4,4,2) + 0.1 # c(bottom, left, top, right)
par(mar = mar.default + c(0, 1, 0, 0)) 
title=paste(mainTitle, " (", clusterNum, " clusters)", sep="")
plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
     main = title, font.main=2, cex.lab=1.6, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
       lty=1, lwd=3, cex=1.3, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
       seg.len = 0.3)
digit=ceiling(-log10(p_value)+2)
if (p_value < 2e-16) {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value < 2e-16"),col="blue",font=2,cex=1.5)
} else {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1.5)  
}
dev.off()

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

x <- surdata
x <- x[,which(colnames(x) %in% c(geneList18, geneList12))]

# change from log count to count
x_count <- 2^x

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
dat <- logCPM

rr <- dat
rr <- rr[which(row.names(rr) %in% gene_set),]
r1 <- rr[,which(colnames(rr) %in% res_high)]
r2 <- cbind(r1, rr[,which(colnames(rr) %in% res_low)])
ncol(r1)
ncol(r2) - ncol(r1)
# > ncol(r1)
# [1] 1177
# > ncol(r2) - ncol(r1)
# [1] 3296

# draw the chart
f <- paste(outDir, "/heatmap_3up.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
data <- as.matrix(r2)
my_group <- c(rep(1, 3))
rowSide <- brewer.pal(3, "Set1")[my_group]
#my_group2 <- c(rep(1, 489), rep(2,390))
my_group2 <- c(rep(1, 1177), rep(2,3296))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------------------------
# 2 down
#------------------------------------
# the cases of 30 genes, including 18 genes and 12 genes, for down genes in 12 genes
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata

exp_data <- surdata[,-c(1,2,3)]
exp_data <- exp_data[,colnames(exp_data) %in% c(geneList18, geneList12)]
exp_data[1:4,1:5]
nrow(exp_data)
ncol(exp_data)
# > exp_data[1:4,1:5]
# ADAM10    AOC3   AP1S1     APP    AQP1
# TCGA.OR.A5J1.01 11.8138  8.0279 10.6627 15.6575 10.6215
# TCGA.OR.A5J2.01 12.6810 10.3151 12.6517 16.5306 11.4720
# TCGA.OR.A5J3.01 10.9106  9.2897 13.5362 16.2096  8.7846
# TCGA.OR.A5J5.01  9.8872  9.0461 11.4528 14.5345  9.1724
# > nrow(exp_data)
# [1] 6727
# > ncol(exp_data)
# [1] 30

x <- exp_data

# change from log count to count
x_count <- 2^x
# > x_count[1:4,1:5]
# ADAM10      AOC3     AP1S1      APP      AQP1
# TCGA.OR.A5J1.01 3600.0469  260.9989  1621.036 51686.50 1575.3973
# TCGA.OR.A5J2.01 6566.9147 1273.9562  6434.891 94668.71 2840.6394
# TCGA.OR.A5J3.01 1924.9430  625.8617 11879.611 75783.70  440.9893
# TCGA.OR.A5J5.01  946.9865  528.6247  2803.085 23731.24  576.9890
# > nrow(x_count)
# [1] 6727
# > ncol(x_count)
# [1] 30

# compute CPM
dat <- t(x_count)
CPM <- cpm(dat, log=FALSE)

exp_data <- CPM
exp_data <- t(exp_data)
mean30 <- colMeans(exp_data)
# > mean30
# ADAM10        AOC3       AP1S1         APP        AQP1        ARF1     ATP6AP1       ATP7A       CASP3        CDK1       COX17          CP 
# 7252.3709   2600.3428   3196.4260  59134.3081  16067.4004  30894.1273  13193.3335   1309.5462   2062.9998   2080.0265   1825.1825  10430.6465 
# CYP1A1         DBH        GPC1       GSK3B         JUN    MAP1LC3A        MAPT      MT-CO1        MT1X        PRNP     S100A12     SLC11A2 
# 266.7921   4534.4520   6674.0926   4320.8433  14387.6737   1209.0188   3583.1980 775925.9396   3771.6661   8363.3308    103.0065   4360.4695 
# SNCA        SORD         SP1     TMPRSS6        XAF1        XIAP 
# 860.8214   7855.0815   7069.0119    702.0836   2211.8583   3753.9499     

result_list <- list()
list_name <- c(1:30)
for (i in 1:30) {
  gene <- colnames(exp_data)[i]
  high <- rownames(exp_data)[which(exp_data[,i] >= mean30[i])]
  low <- rownames(exp_data)[which(exp_data[,i] < mean30[i])]
  result_list[[list_name[i]]] <- list(gene = gene, high = high, low = low)
}

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y <- y[which(y$frequency > 1),]
y_plus <- y
y_plus$UpDown <- "Varied"
y_plus$UpDown <- ifelse(y_plus$numUp > y_plus$numDown, "Up", y_plus$UpDown)
y_plus$UpDown <- ifelse(y_plus$numUp < y_plus$numDown, "Down", y_plus$UpDown)
y_plus$percent <- ifelse(y_plus$UpDown == "Up", y_plus$numUp/y_plus$frequency, y_plus$numDown/y_plus$frequency)
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
y_plus <- y_plus[which(y_plus$gene %in% uni_gene$x),]

# gene_set <- c("ATP7A", "CP", "APP", "MAPT", "DBH", "TMPRSS6", "AOC3", "ADAM10", "SP1", "XIAP")
# gene_set <- c("CDK1", "AP1S1", "CASP3")
# gene_set <- y_plus[which(y_plus$UpDown == "Up" & y_plus$percent >= 1),]$gene
# gene_set <- y_plus[which(y_plus$frequency > 7),]$gene # top 5
# gene_set <- y_plus[which(y_plus$frequency > 7 & y_plus$UpDown == "Up"),]$gene
gene_set <- y_plus[which(y_plus$frequency > 7 & y_plus$UpDown == "Down"),]$gene
gene_set
result_list[[18]]$gene
# > gene_set
# [1] "MAP1LC3A" "SNCA"    
# > result_list[[18]]$gene
# [1] "MAP1LC3A"
# all genes up
res_high <- result_list[[18]]$high
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    res_high <- intersect(res_high, result_list[[i]]$high) 
    print(length(res_high))
  }
}
# [1] 2196
# [1] 562

# genes low
# > result_list[[18]]$gene
# [1] "MAP1LC3A"
res_low <- result_list[[18]]$low
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% gene_set) {
    res_low <- intersect(res_low, result_list[[i]]$low)  
    print(length(res_low))
  }
}
# [1] 4531
# [1] 3993

# tumour classification & survival analysis
surdata$group <- 3
for (i in 1:nrow(surdata)) {
  if(row.names(surdata)[i] %in% res_high) {
    surdata$group[i] <- 1
  }
  if(row.names(surdata)[i] %in% res_low) {
    surdata$group[i] <- 2
  }
}
nrow(surdata[which(surdata$group == 1),]) # 3 genes up
nrow(surdata[which(surdata$group == 2),]) # 3 genes down
# > nrow(surdata[which(surdata$group == 1),]) # 3 genes up
# [1] 562
# > nrow(surdata[which(surdata$group == 2),]) # 3 genes down
# [1] 3993
surdata <- surdata[which(surdata$group %in% c(1,2)),]
clusterNum = length(unique(surdata$group))
time <- surdata$OS.time
status <- surdata$OS
group <- surdata$group
dataset = list(time, status, x = group)  
surv = survfit(Surv(time, status) ~ x, dataset)
mainTitle <- "Survival analysis"
if(clusterNum>1){
  sdf=NULL
  sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
  print(sdf)
  p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
}else{
  cat("There is only one cluster in the group")
  p_value=1
}

# survival chart
# draw the chart
f <- paste(outDir, "/sur_chart_2down.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
# myCol <- wes_palette("Zissou1")
myCol <- wes_palette("Darjeeling2")

# Graph 1
mar.default <- c(5,4,4,2) + 0.1 # c(bottom, left, top, right)
par(mar = mar.default + c(0, 1, 0, 0)) 
title=paste(mainTitle, " (", clusterNum, " clusters)", sep="")
plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
     main = title, font.main=2, cex.lab=1.6, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
       lty=1, lwd=3, cex=1.3, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
       seg.len = 0.3)
digit=ceiling(-log10(p_value)+2)
if (p_value < 2e-16) {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value < 2e-16"),col="blue",font=2,cex=1.5)
} else {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1.5)  
}
dev.off()

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# > surdata[1:5,1:4]
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# > nrow(surdata)
# [1] 6727
# > ncol(surdata)
# [1] 60

x <- surdata
x <- x[,which(colnames(x) %in% c(geneList18, geneList12))]

# change from log count to count
x_count <- 2^x

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
dat <- logCPM

rr <- dat
rr <- rr[which(row.names(rr) %in% gene_set),]
r1 <- rr[,which(colnames(rr) %in% res_high)]
r2 <- cbind(r1, rr[,which(colnames(rr) %in% res_low)])
ncol(r1)
ncol(r2) - ncol(r1)
# > ncol(r1)
# [1] 562
# > ncol(r2) - ncol(r1)
# [1] 3993

# draw the chart
f <- paste(outDir, "/heatmap_2down.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
data <- as.matrix(r2)
my_group <- c(rep(1, 2))
rowSide <- brewer.pal(3, "Set1")[my_group]
#my_group2 <- c(rep(1, 489), rep(2,390))
my_group2 <- c(rep(1, 562), rep(2,3993))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()
#-----------------------------------------

#================================================================

#================================================================
# (28) Process SNV data, evaluate SNVs
# run on server
#================================================================

#-----------------------------------------
# Process SNV data

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]
subtypes
# subtypes
# [1] "ACC"  "AML"  "BLCA" "BRCA" "COAD" "ESCA" "GBM"  "KICH" "KIRC" "KIRP"
# [11] "LGG"  "LIHC" "LUAD" "LUSC" "OVCA" "PAAD" "PCPG" "PRAD" "SKCM" "STAD"
# [21] "THCA" "UCEC" "UCS"

n <- length(subtypes)
outDir <- paste(rootDir, "/Data/Output", sep = "")
for (i in 1:n) {
  cancertype <- subtypes[i]
  
  # load the data
  fileName<-paste0(rootDir, "/Data/SNV/", cancertype, "_maf_data.IdTrans.tsv")
  dat<-as.data.frame(fread(fileName))
  # nrow(dat)
  # dat[1:5,1:10]
  # [1] 11867
  # Center NCBI_Build Chromosome Start_Position End_Position Strand
  # 1      .     GRCh37         10       88419681     88419681      +
  #   2      .     GRCh37         12      133360652    133360652      +
  #   3      .     GRCh37         12        9760409      9760409      +
  #   4      .     GRCh37         14       21991730     21991730      +
  #   5      .     GRCh37         15       40587141     40587141      +
  #   Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1
  # 1      Missense_Mutation          SNP                G                 G
  # 2                 Intron          SNP                C                 C
  # 3      Missense_Mutation          SNP                C                 C
  # 4      Missense_Mutation          SNP                C                 C
  # 5                 Silent          SNP                C                 C
  
  # > unique(dat$Variant_Classification)
  # [1] "Missense_Mutation"      "Intron"                 "Silent"
  # [4] "Frame_Shift_Del"        "3'UTR"                  "Splice_Site"
  # [7] "Nonsense_Mutation"      "RNA"                    "5'UTR"
  # [10] "Frame_Shift_Ins"        "In_Frame_Del"           "3'Flank"
  # [13] "In_Frame_Ins"           "5'Flank"                "Nonstop_Mutation"
  # [16] "Translation_Start_Site"
  
  # Only get seven types of mutation: Missense_Mutation, Nonsense_Mutation, Frame_Shift_Ins, 
  # Splice_Site, Frame_Shift_Del, In_Frame_Del, and In_Frame_Ins. We call these types of mutations deleterious mutations.
  dat <- dat[which(dat$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", 
                                                     "Splice_Site", "Frame_Shift_Del", "In_Frame_Del", "In_Frame_Ins")),]
  # > nrow(dat)
  # unique(dat$Variant_Classification)
  # [1] 7682
  # [1] "Missense_Mutation" "Frame_Shift_Del"   "Splice_Site"
  # [4] "Nonsense_Mutation" "Frame_Shift_Ins"   "In_Frame_Del"
  # [7] "In_Frame_Ins"
  
  # > dat[1:5,38:43]
  # t_depth t_ref_count t_alt_count n_depth n_ref_count n_alt_count
  # 1     143         133          10     140         140           0
  # 3     279         166         113     187         186           0
  # 4     101          54          47      79          78           1
  # 6     191         121          70     113         113           0
  # 7     168         122          45     113         112           0
  
  dat$VAF <- dat$t_alt_count/dat$t_depth
  # dat[1:5,c(38:43,118)]
  # t_depth t_ref_count t_alt_count n_depth n_ref_count n_alt_count        VAF
  # 1     143         133          10     140         140           0 0.06993007
  # 3     279         166         113     187         186           0 0.40501792
  # 4     101          54          47      79          78           1 0.46534653
  # 6     191         121          70     113         113           0 0.36649215
  # 7     168         122          45     113         112           0 0.26785714
  
  dat <- dat[which(dat$FILTER == "PASS"),]
  
  fwrite(dat, paste0(outDir, "/SNV/", cancertype, "_maf_data.IdTrans.deleterious.csv"))
}

#-----------------------------------------
# Evaluate SNVs
# Remember to remove SNVs with VAF < 0.05

threshold <- 0.05

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get 30 genes

outDir <- paste(rootDir, "/Data/Output", sep = "")

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# nrow(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# [1] 35

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]

geneList30 <- c(geneList18, geneList12)
geneList30
# [1] "CDK1"     "AP1S1"    "CASP3"    "TMPRSS6"  "GSK3B"    "APP"
# [7] "COX17"    "XIAP"     "ARF1"     "GPC1"     "SORD"     "ATP7A"
# [13] "SP1"      "MT-CO1"   "SLC11A2"  "ATP6AP1"  "ADAM10"   "CP"
# [19] "MAP1LC3A" "SNCA"     "MAPT"     "JUN"      "CYP1A1"   "AOC3"
# [25] "PRNP"     "DBH"      "S100A12"  "AQP1"     "MT1X"     "XAF1"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get sample size of cancer types

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

n <- length(subtypes)
r <- matrix(NA, nrow = n, ncol = 2)
colnames(r) <- c("cancer_type", "sample_size")
for (i in 1:n) {
  cancertype <- subtypes[i]
  
  # load the data
  # for total number of samples
  fileName<-paste0(rootDir, "/Data/SNV/", cancertype, "_maf_data.IdTrans.tsv")
  dat<-as.data.frame(fread(fileName))
  
  r[i,1] <- cancertype
  r[i,2] <- length(unique(dat$barcode16))
  # r[i,2] <- length(unique(dat$Tumor_Sample_Barcode))

}

r <- as.data.frame(r)
r[,2] <- as.numeric(r[,2])
r
# cancer_type sample_size
# 1          ACC          92
# 2          AML          85
# 3         BLCA         411
# 4         BRCA        1026
# 5         COAD         407
# 6         ESCA         185
# 7          GBM         403
# 8         KICH          66
# 9         KIRC         370
# 10        KIRP         282
# 11         LGG         526
# 12        LIHC         365
# 13        LUAD         567
# 14        LUSC         485
# 15        OVCA         412
# 16        PAAD         178
# 17        PCPG         184
# 18        PRAD         498
# 19        SKCM         468
# 20        STAD         439
# 21        THCA         500
# 22        UCEC         531
# 23         UCS          57

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# draw heat map

# geneList30

# r
r$type <- paste0(r$cancer_type, " (n=", r$sample_size, ")")

nCancers <- nrow(r)
nGenes <- length(geneList30)
dat <- matrix(0, nrow = nGenes, ncol = nCancers)
colnames(dat) <- r$type
row.names(dat) <- geneList30
for (i in 1:nCancers) {
  cancertype <- r$cancer_type[i]
  
  # load the data
  fileName<-paste0(outDir, "/SNV/", cancertype, "_maf_data.IdTrans.deleterious.csv")
  mut<-as.data.frame(fread(fileName))
  
  mut <- mut[which(mut$VAF >= threshold),]
  
  for (j in 1:nGenes) {
    gene <- geneList30[j]
    dat[j, i] <- length(unique(mut[which(mut$Hugo_Symbol == gene),"barcode16"]))
  }
}

dat2 <- melt(dat)
dat2$percent <- dat2$value/parse_number(as.character(dat2$Var2)) 
colnames(dat2) <- c("Gene", "Cancer_Type", "value", "Mutation_Frequency")
dat2[which(dat2$value==0),"value"] <- ""

# draw the chart
f <- paste(outDir, "/SNV_heatmap.pdf", sep = "")
pdf(file = f, width = 7, height =  7)
ggplot(dat2, aes(Cancer_Type, Gene)) +
  geom_tile(aes(fill = Mutation_Frequency)) +
  geom_text(aes(label = value)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get mutation data

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

SNV_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del",
               "In_Frame_Ins")

# geneList30

# r - list of cancer types

# get all mutation info
n <- length(subtypes)
for (i in 1:n) {
  cancertype <- subtypes[i]
  # Read data
  dat <- as.data.frame(fread(paste0(outDir, "/SNV/", cancertype, "_maf_data.IdTrans.deleterious.csv")))
  dat <- dat[which(dat$VAF >= threshold),]
  # Only get necessary data
  dat <- dat[,c("Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "Allele", "barcode16", "Hugo_Symbol")]
  if(nrow(dat) > 0) {
    dat$CancerType <- cancertype
    if(i == 1) {
      res <- dat
    } else {
      res <- rbind(res,dat)
    }  
  }
}

# Mutations for all genes
saveRDS(res, paste0(outDir, "/SNV/SNVAllGenes_005.rds"))
# res <- readRDS(paste0(outDir, "/SNV/SNVAllGenes_005.rds"))

# 30 genes
saveRDS(res[which(res$Hugo_Symbol %in% geneList30),], 
        paste0(outDir, "/SNV/SNV30Genes_005.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Draw bar chart for Variant_Classification

# data of 30 genes
mut_30 <- readRDS(paste0(outDir, "/SNV/SNV30Genes_005.rds"))
# > head(mut_30)
# Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele2
# 374          Frame_Shift_Del          DEL                C                 -
#   680        Missense_Mutation          SNP                G                 T
# 5808         Frame_Shift_Del          DEL                G                 -
#   5963       Missense_Mutation          SNP                C                 G
# 6505       Missense_Mutation          SNP                G                 C
# 12930      Missense_Mutation          SNP                C                 G
# Allele        barcode16 Hugo_Symbol CancerType
# 374        - TCGA-OR-A5J5-01A     TMPRSS6        ACC
# 680        T TCGA-OR-A5J8-01A         DBH        ACC
# 5808       - TCGA-OR-A5LJ-01A       AP1S1        ACC
# 5963       G TCGA-OR-A5LL-01A          CP        ACC
# 6505       C TCGA-P6-A5OH-01A       AP1S1        ACC
# 12930      G TCGA-2F-A9KO-01A      CYP1A1       BLCA
# > table(mut_30$Variant_Classification)
# 
# Frame_Shift_Del   Frame_Shift_Ins      In_Frame_Del      In_Frame_Ins
# 113                23                11                 1
# Missense_Mutation Nonsense_Mutation       Splice_Site
# 1625               105                38

variant_data <- as.data.frame(table(mut_30$Variant_Classification))
variant_data <- variant_data[order(variant_data$Freq, decreasing = FALSE),]
# > variant_data
# Var1 Freq
# 4      In_Frame_Ins    1
# 3      In_Frame_Del   11
# 2   Frame_Shift_Ins   23
# 7       Splice_Site   38
# 6 Nonsense_Mutation  105
# 1   Frame_Shift_Del  113
# 5 Missense_Mutation 1625

# Specific color for each bar? Use a well known palette
f <- paste(outDir, "/SNV_variant_classification.pdf", sep = "")
pdf(file = f, width = 4, height =  4)
# coul <- brewer.pal(7, "Set2") 
# > coul
# [1] "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
# adjust to the maximum of either the default 
# or a figure based on the maximum length
# coul <- c("brown", "yellow", "violet", "orange", "red", "blue", "green4")
# mutType <- c("In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Missense_Mutation")
coul <- c("brown", "yellow", "violet", "orange", "red", "blue", "green4")
ylabels <- c("Missense_Mutation")
par(mar=c(5.1, max(4.1,max(nchar(ylabels))/1.8) ,4.1 ,2.1))
for (i in 1:nrow(variant_data)) {
  if(variant_data$Var1[i] == "In_Frame_Ins") {
    variant_data$Color[i] <- "brown"
  } else if(variant_data$Var1[i] == "In_Frame_Del") {
    variant_data$Color[i] <- "yellow"
  } else if(variant_data$Var1[i] == "Frame_Shift_Ins") {
    variant_data$Color[i] <- "violet"
  } else if(variant_data$Var1[i] == "Splice_Site") {
    variant_data$Color[i] <- "orange"
  } else if(variant_data$Var1[i] == "Nonsense_Mutation") {
    variant_data$Color[i] <- "red"
  } else if(variant_data$Var1[i] == "Frame_Shift_Del") {
    variant_data$Color[i] <- "blue"
  } else if(variant_data$Var1[i] == "Missense_Mutation") {
    variant_data$Color[i] <- "green4"
  }
}
barplot(height=variant_data$Freq, names=variant_data$Var1, col=variant_data$Color, horiz=T, las=1, border=NA)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Draw bar chart for top 10 genes

# data of 30 genes
mut_30 <- readRDS(paste0(outDir, "/SNV/SNV30Genes_005.rds"))
# > head(mut_30)
# Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele2
# 374          Frame_Shift_Del          DEL                C                 -
#   680        Missense_Mutation          SNP                G                 T
# 5808         Frame_Shift_Del          DEL                G                 -
#   5963       Missense_Mutation          SNP                C                 G
# 6505       Missense_Mutation          SNP                G                 C
# 12930      Missense_Mutation          SNP                C                 G
# Allele        barcode16 Hugo_Symbol CancerType
# 374        - TCGA-OR-A5J5-01A     TMPRSS6        ACC
# 680        T TCGA-OR-A5J8-01A         DBH        ACC
# 5808       - TCGA-OR-A5LJ-01A       AP1S1        ACC
# 5963       G TCGA-OR-A5LL-01A          CP        ACC
# 6505       C TCGA-P6-A5OH-01A       AP1S1        ACC
# 12930      G TCGA-2F-A9KO-01A      CYP1A1       BLCA

gene_frq <- as.data.frame(table(mut_30$Hugo_Symbol))
gene_frq <- gene_frq[order(gene_frq$Freq, decreasing = TRUE),]
total_mut <- nrow(mut_30)
# > total_mut
# [1] 1916
gene_frq$percent <- gene_frq$Freq/total_mut
gene_frq$percent <- paste0(round(gene_frq$percent*100,0), "%")
gene_frq <- gene_frq[1:10,]
# Var1 Freq percent
# 8    ATP7A  223     12%
# 12      CP  159      8%
# 19    MAPT  123      6%
# 4      APP  122      6%
# 27 TMPRSS6  112      6%
# 14     DBH  108      6%
# 1   ADAM10   98      5%
# 2     AOC3   89      5%
# 26     SP1   87      5%
# 13  CYP1A1   75      4%

SNV_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del",
               "In_Frame_Ins")

dat <- gene_frq
# total_mut <- sum(dat$Freq)
# dat$percent <- paste0(round((dat$Freq/total_mut)*100, 1), "%")
dat$Missense_Mutation <- 0
dat$Nonsense_Mutation <- 0
dat$Frame_Shift_Del <- 0
dat$Splice_Site <- 0
dat$Frame_Shift_Ins <- 0
dat$In_Frame_Del <- 0
dat$In_Frame_Ins <- 0
for (i in 1:10) {
  mut <- mut_30[which(mut_30$Hugo_Symbol %in% dat$Var1[i]),]
  n <- nrow(mut)
  for (j in 1:n) {
    if(mut$Variant_Classification[j] == "Missense_Mutation") {
      dat$Missense_Mutation[i] <- dat$Missense_Mutation[i] + 1
      
    } else if(mut$Variant_Classification[j] == "Nonsense_Mutation") {
      dat$Nonsense_Mutation[i] <- dat$Nonsense_Mutation[i] + 1
      
    } else if(mut$Variant_Classification[j] == "Frame_Shift_Del") {
      dat$Frame_Shift_Del[i] <- dat$Frame_Shift_Del[i] + 1
      
    } else if(mut$Variant_Classification[j] == "Splice_Site") {
      dat$Splice_Site[i] <- dat$Splice_Site[i] + 1
      
    } else if(mut$Variant_Classification[j] == "Frame_Shift_Ins") {
      dat$Frame_Shift_Ins[i] <- dat$Frame_Shift_Ins[i] + 1
      
    } else if(mut$Variant_Classification[j] == "In_Frame_Del") {
      dat$In_Frame_Del[i] <- dat$In_Frame_Del[i] + 1
      
    } else if(mut$Variant_Classification[j] == "In_Frame_Ins") {
      dat$In_Frame_Ins[i] <- dat$In_Frame_Ins[i] + 1
    }
    
  }
}
dat
# dat
# Var1 Freq percent Missense_Mutation Nonsense_Mutation Frame_Shift_Del
# 8    ATP7A  223     12%               199                14               5
# 12      CP  159      8%               126                 8              18
# 19    MAPT  123      6%               105                 3              11
# 4      APP  122      6%               105                 6               5
# 27 TMPRSS6  112      6%               100                 4               6
# 14     DBH  108      6%                97                 6               3
# 1   ADAM10   98      5%                81                 4               8
# 2     AOC3   89      5%                78                 2               6
# 26     SP1   87      5%                73                 9               2
# 13  CYP1A1   75      4%                66                 3               2
# Splice_Site Frame_Shift_Ins In_Frame_Del In_Frame_Ins
# 8            3               2            0            0
# 12           5               1            1            0
# 19           3               1            0            0
# 4            4               1            1            0
# 27           1               0            1            0
# 14           2               0            0            0
# 1            1               2            2            0
# 2            0               3            0            0
# 26           1               1            1            0
# 13           1               2            0            1

dat2 <- data.frame(Gene = rep(dat$Var1, 7), TotalNumMut = rep(dat$Freq, 7), Per = rep(dat$percent, 7), 
                   MutType = c(rep(SNV_types[1], 10),
                               rep(SNV_types[2], 10),
                               rep(SNV_types[3], 10),
                               rep(SNV_types[4], 10),
                               rep(SNV_types[5], 10),
                               rep(SNV_types[6], 10),
                               rep(SNV_types[7], 10)),
                   NumMut = c(dat$Missense_Mutation,
                              dat$Nonsense_Mutation,
                              dat$Frame_Shift_Del,
                              dat$Splice_Site,
                              dat$Frame_Shift_Ins,
                              dat$In_Frame_Del,
                              dat$In_Frame_Ins),
                   Color = "")
# for (i in 1:nrow(dat2)) {
#   if(dat2$MutType[i] == "In_Frame_Ins") {
#     dat2$Color[i] <- "#E5C494"
#   } else if(dat2$MutType[i] == "In_Frame_Del") {
#     dat2$Color[i] <- "#FFD92F"
#   } else if(dat2$MutType[i] == "Frame_Shift_Ins") {
#     dat2$Color[i] <- "#A6D854"
#   } else if(dat2$MutType[i] == "Splice_Site") {
#     dat2$Color[i] <- "#E78AC3"
#   } else if(dat2$MutType[i] == "Nonsense_Mutation") {
#     dat2$Color[i] <- "#8DA0CB"
#   } else if(dat2$MutType[i] == "Frame_Shift_Del") {
#     dat2$Color[i] <- "#FC8D62"
#   } else if(dat2$MutType[i] == "Missense_Mutation") {
#     dat2$Color[i] <- "#66C2A5"
#   }
# }
# SNV_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
#                "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del",
#                "In_Frame_Ins")
# coul <- c("brown", "yellow", "violet", "orange", "red", "blue", "green4")
# mutType <- c("In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Missense_Mutation")
for (i in 1:nrow(dat2)) {
  if(dat2$MutType[i] == "In_Frame_Ins") {
    dat2$Color[i] <- "brown"
  } else if(dat2$MutType[i] == "In_Frame_Del") {
    dat2$Color[i] <- "yellow"
  } else if(dat2$MutType[i] == "Frame_Shift_Ins") {
    dat2$Color[i] <- "violet"
  } else if(dat2$MutType[i] == "Splice_Site") {
    dat2$Color[i] <- "orange"
  } else if(dat2$MutType[i] == "Nonsense_Mutation") {
    dat2$Color[i] <- "red"
  } else if(dat2$MutType[i] == "Frame_Shift_Del") {
    dat2$Color[i] <- "blue"
  } else if(dat2$MutType[i] == "Missense_Mutation") {
    dat2$Color[i] <- "green4"
  }
}
save(dat2, file=paste0(outDir, 'dat2.RData'))
# > head(dat2)
# Gene TotalNumMut Per           MutType NumMut Color
# 1   ATP7A         223 12% Missense_Mutation    199 green
# 2      CP         159  8% Missense_Mutation    126 green
# 3    MAPT         123  6% Missense_Mutation    105 green
# 4     APP         122  6% Missense_Mutation    105 green
# 5 TMPRSS6         112  6% Missense_Mutation    100 green
# 6     DBH         108  6% Missense_Mutation     97 green

# chart
f <- paste(outDir, "/SNV_Top10Genes.pdf", sep = "")
pdf(file = f, width = 4, height =  4)
## set the levels in order we want

load(paste0(outDir, 'dat2.RData'))
# load("R:/TTRP/Vu/002NetworkAnalysis/Data/Output/dat2.RData")

dat3 <- within(dat2, Gene <- factor(Gene, levels=c("CYP1A1", "SP1", "AOC3",
                                                   "ADAM10", "DBH", "TMPRSS6",
                                                   "APP", "MAPT", "CP", "ATP7A")))
# dat3 <- within(dat3, MutType <- factor(MutType, levels=c("In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins",
#                                                          "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Missense_Mutation")))
dat3 <- within(dat3, MutType <- factor(MutType, levels=c("In_Frame_Ins", "Frame_Shift_Ins", "Splice_Site",
                                                         "Nonsense_Mutation", "In_Frame_Del", "Missense_Mutation",
                                                         "Frame_Shift_Del")))
# par(mar=c(5.1, 4.1 ,4.1 ,4.1)) # bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1)
# coul <- c("brown", "yellow", "violet", "orange", "red", "blue", "green4")
# mutType <- c("In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Missense_Mutation")
lim <- c(0, sum(dat3[which(dat3$Gene == "ATP7A"),"NumMut"]) + 50)
ggplot(aes(y=NumMut, x=Gene, fill = MutType), data = dat3) +
  # geom_bar(stat = 'identity', position = "stack") +
  geom_bar(stat = 'identity') +
  scale_y_continuous(expand = c(0, 0), limits = lim) + 
  geom_text(data = dat3[1:10,],
            hjust = -0.5,
            aes(x= Gene, y = TotalNumMut, label = Per),
            inherit.aes = FALSE) +
  coord_flip() +
  scale_fill_manual('Variant Classification', values = c("In_Frame_Ins" = "brown", "In_Frame_Del" = "yellow", "Frame_Shift_Ins" = "violet",
                                                         "Splice_Site" = "orange", "Nonsense_Mutation" = "red", "Frame_Shift_Del" = "blue",
                                                         "Missense_Mutation" = "green4")) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(), legend.position="none")
  # theme(plot.margin = margin(1,1,1.5,1.2, "cm"))
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# water fall chart

# Load the data and create maf file
subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

SNV_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del",
               "In_Frame_Ins")

# 30 genes
# geneList30

# get all mutation info
n <- length(subtypes)
outDir <- paste(rootDir, "/Data/Output", sep = "")
threshold <- 0.05

for (i in 1:n) {
  cancertype <- subtypes[i]
  # Read data
  dat <- as.data.frame(fread(paste0(outDir, "/SNV/", cancertype, "_maf_data.IdTrans.deleterious.csv")))
  dat <- dat[which(dat$VAF >= threshold),]
  
  # Get 30 genes
  dat <- dat[which(dat$Hugo_Symbol %in% geneList30),]
  
  # Only get necessary data
  
  # Required columns: Hugo_Symbol	Entrez_Gene_Id	Center	NCBI_Build	Chromosome
  # Start_Position	End_position	Strand	Variant_Classification	Variant_Type
  # Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	Tumor_Sample_Barcode	Protein_Change
  # i_TumorVAF_WU	i_transcript_name
  
  # Note
  # VAF, need to *100
  # Use Transcript_ID instead of i_transcript_name
  dat_mut <- dat[,c("Hugo_Symbol",	"entrez",	"Center",	"NCBI_Build",	"Chromosome",
    "Start_Position",	"End_Position",	"Strand",	"Variant_Classification",	"Variant_Type",
    "Reference_Allele",	"Tumor_Seq_Allele1",	"Tumor_Seq_Allele2",	"Tumor_Sample_Barcode",	"HGVSp_Short", 
    "VAF",	"Transcript_ID")]
  colnames(dat_mut) <- c("Hugo_Symbol",	"Entrez_Gene_Id",	"Center",	"NCBI_Build",	"Chromosome",
                     "Start_Position",	"End_position",	"Strand",	"Variant_Classification",	"Variant_Type",
                     "Reference_Allele",	"Tumor_Seq_Allele1",	"Tumor_Seq_Allele2",	"Tumor_Sample_Barcode",	"Protein_Change",
                     "i_TumorVAF_WU",	"i_transcript_name")
  dat_mut$i_TumorVAF_WU <- dat_mut$i_TumorVAF_WU*100
  
  dat_clinical <- dat[,c("Tumor_Sample_Barcode", "cancer_types")]
  colnames(dat_clinical) <- c("Tumor_Sample_Barcode", "Cancer_Type")
  
  if(nrow(dat_mut) > 0) {
    dat_clinical$Cancer_Type <- cancertype
    if(i == 1) {
      res_mut <- dat_mut
      res_clinical <- dat_clinical
    } else {
      res_mut <- rbind(res_mut,dat_mut)
      res_clinical <- rbind(res_clinical, dat_clinical)
    }  
  }
}

write.table(res_mut, paste0(outDir, "/pan_cancer_30genes.maf"), quote = FALSE, row.names = F, sep="\t")
write.table(res_clinical, paste0(outDir, "/pan_cancer_annot_30genes.tsv"), quote = FALSE, row.names = F, sep="\t")

# draw chart
# f <- paste(outDir, "/SNV_Waterfall.pdf", sep = "")
# pdf(file = f, width = 12, height =  12)
f <- paste(outDir, "/SNV_Waterfall.png", sep = "")
png(file = f, width = 900, height =  600)
pan.maf <- paste0(outDir, "/pan_cancer_30genes.maf")
pan.clin <- paste0(outDir, "/pan_cancer_annot_30genes.tsv")
pan <- read.maf(maf = pan.maf, clinicalData = pan.clin)
#Basic onocplot
# oncoplot(maf = pan, top = 3)
#Changing colors for variant classifications (You can use any colors, here in this example we will use a color palette from RColorBrewer)
# for (i in 1:nrow(dat2)) {
#   if(dat2$MutType[i] == "In_Frame_Ins") {
#     dat2$Color[i] <- "#E5C494"
#   } else if(dat2$MutType[i] == "In_Frame_Del") {
#     dat2$Color[i] <- "#FFD92F"
#   } else if(dat2$MutType[i] == "Frame_Shift_Ins") {
#     dat2$Color[i] <- "#A6D854"
#   } else if(dat2$MutType[i] == "Splice_Site") {
#     dat2$Color[i] <- "#E78AC3"
#   } else if(dat2$MutType[i] == "Nonsense_Mutation") {
#     dat2$Color[i] <- "#8DA0CB"
#   } else if(dat2$MutType[i] == "Frame_Shift_Del") {
#     dat2$Color[i] <- "#FC8D62"
#   } else if(dat2$MutType[i] == "Missense_Mutation") {
#     dat2$Color[i] <- "#66C2A5"
#   }
# }
# col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
# names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
#                'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
# coul <- c("brown", "yellow", "violet", "orange", "red", "blue", "green4")
# mutType <- c("In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Missense_Mutation")
col = c("brown", "yellow", "violet", "orange", "red", "blue", "green4", "black")
names(col) = c("In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Missense_Mutation", "Multi-Hit")
#Color coding for classification; try getAnnotations(x = pan) to see available annotations.
# fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
mut_30 <- readRDS(paste0(outDir, "/SNV/SNV30Genes_005.rds"))
length(unique(mut_30$CancerType))
# length(unique(mut_30$CancerType))
# [1] 22
# > unique(mut_30$CancerType)
# [1] "ACC"  "BLCA" "BRCA" "COAD" "ESCA" "GBM"  "KICH" "KIRC" "KIRP" "LGG"
# [11] "LIHC" "LUAD" "LUSC" "OVCA" "PAAD" "PCPG" "PRAD" "SKCM" "STAD" "THCA"
# [21] "UCEC" "UCS"
fabcolors = c(RColorBrewer::brewer.pal(n = 12, name = 'Paired'),RColorBrewer::brewer.pal(n = 10, name = 'Set3'))
names(fabcolors) = unique(mut_30$CancerType)
fabcolors = list(Cancer_Type = fabcolors)
oncoplot(maf = pan, colors = col, clinicalFeatures = 'Cancer_Type', annotationColor = fabcolors, top = 10, sepwd_samples = 0,
         sortByAnnotation = TRUE,
         groupAnnotationBySize = FALSE,
         annotationOrder = NULL,
         sortByMutation = FALSE)
dev.off()

# Check percentages
# data of 30 genes
mut_30 <- readRDS(paste0(outDir, "/SNV/SNV30Genes_005.rds"))

gene_frq <- as.data.frame(table(mut_30$Hugo_Symbol))
gene_frq <- gene_frq[order(gene_frq$Freq, decreasing = TRUE),]
gene_frq <- gene_frq[1:10,]
# > gene_frq
# Var1 Freq
# 8    ATP7A  223
# 12      CP  159
# 19    MAPT  123
# 4      APP  122
# 27 TMPRSS6  112
# 14     DBH  108
# 1   ADAM10   98
# 2     AOC3   89
# 26     SP1   87
# 13  CYP1A1   75

mut_10 <- mut_30[which(mut_30$Hugo_Symbol %in% gene_frq$Var1),]
gene_frq$per <- 0
for (i in 1:10) {
  mut <- mut_10[which(mut_10$Hugo_Symbol == gene_frq$Var1[i]),]
  gene_frq$per[i] <- round(length(unique(mut$barcode16))/length(unique(mut_30$barcode16))*100,0)
}
gene_frq
# gene_frq
# Var1 Freq per
# 8    ATP7A  223  17
# 12      CP  159  13
# 19    MAPT  123  11
# 4      APP  122  11
# 27 TMPRSS6  112  10
# 14     DBH  108  10
# 1   ADAM10   98   9
# 2     AOC3   89   8
# 26     SP1   87   8
# 13  CYP1A1   75   7

#-----------------------------------------

#================================================================

#================================================================
# (29) Check location of SNV, gene ATP7A
# run on server
#================================================================
outDir <- paste0(rootDir, "/Data/Output")
threshold <- 0.05
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get mutation data with location

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

SNV_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del",
               "In_Frame_Ins")

# gene ATP7A
gene <- "ATP7A"

# get all mutation info
n <- length(subtypes)
for (i in 1:n) {
  cancertype <- subtypes[i]
  # Read data
  dat <- as.data.frame(fread(paste0(outDir, "/SNV/", cancertype, "_maf_data.IdTrans.deleterious.csv")))
  dat <- dat[which(dat$VAF >= threshold),]
  # Only get necessary data
  dat <- dat[,c("Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "Allele", "barcode16", "Hugo_Symbol",
                "Chromosome", "Start_Position", "End_Position", "Exon_Number", "all_effects", "EXON", "INTRON", "IMPACT")]
  if(nrow(dat) > 0) {
    dat$CancerType <- cancertype
    if(i == 1) {
      res <- dat
    } else {
      res <- rbind(res,dat)
    }  
  }
}

# ATP7A
saveRDS(res[which(res$Hugo_Symbol == gene),], 
        paste0(outDir, "/SNV/SNVATP7A_005.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Draw bar chart for Variant_Classification

# data of ATP7A
mut_ATP7A <- readRDS(paste0(outDir, "/SNV/SNVATP7A_005.rds"))
# > table(mut_ATP7A$Variant_Classification)
# 
# Frame_Shift_Del   Frame_Shift_Ins Missense_Mutation Nonsense_Mutation
# 5                 2               199                14
# Splice_Site
# 3
# > nrow(mut_ATP7A)
# [1] 223

variant_data <- as.data.frame(table(mut_ATP7A$Variant_Classification))
variant_data <- variant_data[order(variant_data$Freq, decreasing = FALSE),]

# Specific color for each bar? Use a well known palette
f <- paste(outDir, "/SNV_variant_classification_ATP7A.pdf", sep = "")
pdf(file = f, width = 4, height =  4)
coul <- c("brown", "yellow", "violet", "orange", "red", "blue", "green4")
ylabels <- c("Missense_Mutation")
par(mar=c(5.1, max(4.1,max(nchar(ylabels))/1.8) ,4.1 ,2.1))
for (i in 1:nrow(variant_data)) {
  if(variant_data$Var1[i] == "In_Frame_Ins") {
    variant_data$Color[i] <- "brown"
  } else if(variant_data$Var1[i] == "In_Frame_Del") {
    variant_data$Color[i] <- "yellow"
  } else if(variant_data$Var1[i] == "Frame_Shift_Ins") {
    variant_data$Color[i] <- "violet"
  } else if(variant_data$Var1[i] == "Splice_Site") {
    variant_data$Color[i] <- "orange"
  } else if(variant_data$Var1[i] == "Nonsense_Mutation") {
    variant_data$Color[i] <- "red"
  } else if(variant_data$Var1[i] == "Frame_Shift_Del") {
    variant_data$Color[i] <- "blue"
  } else if(variant_data$Var1[i] == "Missense_Mutation") {
    variant_data$Color[i] <- "green4"
  }
}
barplot(height=variant_data$Freq, names=variant_data$Var1, col=variant_data$Color, horiz=T, las=1, border=NA)
dev.off()

# Location of mutations

# > min(mut_ATP7A$Start_Position)
# [1] 77227117
# > max(mut_ATP7A$End_Position)
# [1] 77302020
# > unique(mut_ATP7A$Chromosome)
# [1] "X"
# > table(mut_ATP7A$Exon_Number)
# 
# . 10/23 11/23 12/23 13/23 14/23 15/23 16/23 17/23 18/23 19/23 20/23 21/23
# 3     9     4     6     9     3    14     6    14     8     6     4     9
# 22/23 23/23  3/23  4/23  5/23  6/23  7/23  8/23  9/23
# 4    15    20    37    17     7    15     3    10

variant_data <- as.data.frame(table(mut_ATP7A$Exon_Number))
variant_data <- variant_data[order(variant_data$Freq, decreasing = FALSE),]

# Specific color for each bar? Use a well known palette
f <- paste(outDir, "/SNV_exon_classification_ATP7A.pdf", sep = "")
pdf(file = f, width = 5, height =  10)
coul <- rep("grey", nrow(variant_data))
# ylabels <- c("Missense_Mutation")
# par(mar=c(5.1, max(4.1,max(nchar(ylabels))/1.8) ,4.1 ,2.1))
# for (i in 1:nrow(variant_data)) {
#   if(variant_data$Var1[i] == "In_Frame_Ins") {
#     variant_data$Color[i] <- "brown"
#   } else if(variant_data$Var1[i] == "In_Frame_Del") {
#     variant_data$Color[i] <- "yellow"
#   } else if(variant_data$Var1[i] == "Frame_Shift_Ins") {
#     variant_data$Color[i] <- "violet"
#   } else if(variant_data$Var1[i] == "Splice_Site") {
#     variant_data$Color[i] <- "orange"
#   } else if(variant_data$Var1[i] == "Nonsense_Mutation") {
#     variant_data$Color[i] <- "red"
#   } else if(variant_data$Var1[i] == "Frame_Shift_Del") {
#     variant_data$Color[i] <- "blue"
#   } else if(variant_data$Var1[i] == "Missense_Mutation") {
#     variant_data$Color[i] <- "green4"
#   }
# }
barplot(height=variant_data$Freq, names=variant_data$Var1, col=variant_data$Color, horiz=T, las=1, border=NA,
        xlab = "Mutation count", ylab = "Exon number")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#================================================================
#================================================================
# (30) CNV & methylation
# run on server
#================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get 30 genes

outDir <- paste(rootDir, "/Data/Output", sep = "")

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]

geneList30 <- c(geneList18, geneList12)
paste(geneList30, collapse = ' ')
# paste(geneList30, collapse = ' ')
# [1] "CDK1 AP1S1 CASP3 TMPRSS6 GSK3B APP COX17 XIAP ARF1 GPC1 SORD ATP7A SP1 MT-CO1 SLC11A2 ATP6AP1 ADAM10 CP MAP1LC3A SNCA MAPT JUN CYP1A1 AOC3 PRNP DBH S100A12 AQP1 MT1X XAF1"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#================================================================

#================================================================
# (31) LGG
# Compare gene expression, copper genes, up genes and down genes, using count and cpm
# heatmap for clustered groups
# identify subset of genes for tumour classification & survival analysis
# run on server
#================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LGG

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
cancertype <- "LGG"
dat <- surdata[which(surdata$cancertype == cancertype),]
nrow(dat)
dat[1:3,1:4]
# [1] 348
# cancertype OS OS.time  AANAT
# TCGA.CS.4938.01        LGG  0    3574 3.3219
# TCGA.CS.4941.01        LGG  1     234 4.7549
# TCGA.CS.4942.01        LGG  1    1335 3.5850

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y[1:4,1:5]
# gene frequency ACC  AML BLCA
# 1  CDK1        18  UP   UP   UP
# 2   ALB        13   0 DOWN    0
# 3 AP1S1        11   0 DOWN    0
# 4 CASP3         9   0    0    0
# Critical copper genes
criticalCopperGenes <- y
upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'UP'),]$gene
downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'DOWN'),]$gene
geneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] %in% c("UP", "DOWN")),]$gene
upGeneList
downGeneList
geneList
# [1] "CDK1"  "CASP3" "XIAP"  "TP53"  "F5"    "ATP7A" "SP1"   "CP"
# [1] "ALB"      "MAP1LC3A" "SNCA"     "CYP1A1"   "MT-CO1"
# [1] "CDK1"     "ALB"      "CASP3"    "MAP1LC3A" "SNCA"     "CYP1A1"
# [7] "XIAP"     "TP53"     "F5"       "ATP7A"    "SP1"      "MT-CO1"
# [13] "CP"

#------------------------------------
# For 3 groups
# select a set of genes which have the same pattern
# e.g., 5 genes while 3 genes high and 2 genes low in a set of samples,
# and 3 genes low and 2 genes high in another set
#------------------------------------
# For 13 genes, identify 2 groups up & down, LGG
# Top5, including both up and down

# Get genes with high node degree
critical_nodes <- read.csv(paste(rootDir, "/Data/LGG/critical_nodes.csv", sep = ""))
copper.critical<- critical_nodes[critical_nodes$Name %in% geneList, ]
copper.critical <- copper.critical[order(copper.critical$K, decreasing = TRUE),]
top5 <- copper.critical[1:5,]
top5
# Name   K Kin Kout TypeI TypeII
# 143  CDK1 180 148   32     0      1
# 445  TP53 151   5  146     0      1
# 109 CASP3  66  25   41     0      1
# 81    ALB  47  18   29     0      1
# 130  SNCA  34  15   19     0      1
# up: CDK1, TP53, CASP3, 
# down: ALB, SNCA

surdata <- dat

# for all 13 genes
exp_data <- surdata[,-c(1,2,3)]
exp_data <- exp_data[,colnames(exp_data) %in% geneList]
exp_data[1:4,1:5]
nrow(exp_data)
ncol(exp_data)
# ALB  ATP7A   CASP3    CDK1      CP
# TCGA.CS.4938.01 5.4263 8.4838 10.3608  6.3923 11.5159
# TCGA.CS.4941.01 5.7808 9.1396 10.8556  9.0084 10.2872
# TCGA.CS.4942.01 4.0875 9.3772 10.6027  8.9425 10.1364
# TCGA.CS.4943.01 3.8074 9.5018 10.4798 12.3250  8.8437
# [1] 348
# [1] 13

x <- exp_data

# change from log count to count
x_count <- 2^x

# compute CPM
dat <- t(x_count)
CPM <- cpm(dat, log=FALSE)

exp_data <- CPM
exp_data <- t(exp_data)
mean13 <- colMeans(exp_data)

result_list <- list()
list_name <- c(1:13)
for (i in 1:13) {
  gene <- colnames(exp_data)[i]
  high <- rownames(exp_data)[which(exp_data[,i] >= mean13[i])]
  low <- rownames(exp_data)[which(exp_data[,i] < mean13[i])]
  result_list[[list_name[i]]] <- list(gene = gene, high = high, low = low)
}

gene_set <- top5$Name
# > gene_set
# [1] "CDK1"  "TP53"  "CASP3" "ALB"   "SNCA"
# > colnames(exp_data)
# [1] "ALB"      "ATP7A"    "CASP3"    "CDK1"     "CP"       "CYP1A1"
# [7] "F5"       "MAP1LC3A" "MT-CO1"   "SNCA"     "SP1"      "TP53"
# [13] "XIAP"
# > result_list[[4]]$gene
# [1] "CDK1"
# up: CDK1, TP53, CASP3, 
# down: ALB, SNCA
# 3 up genes & 2 down genes
res_high <- result_list[[4]]$high # index 4 for CDK1
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% c("CDK1", "TP53", "CASP3")) { # only for top 3
    if(result_list[[i]]$gene %in% c("CDK1", "TP53", "CASP3")) { # 3 up genes
      res_high <- intersect(res_high, result_list[[i]]$high)  
    } else { # 2 down genes
      res_high <- intersect(res_high, result_list[[i]]$low)
    }
    print(length(res_high))
  }
}
# [1] 54
# [1] 54
# [1] 51

# "CDK1", "TP53", "CASP3" low
# > result_list[[4]]$gene
# [1] "CDK1"
res_low <- result_list[[4]]$low # index 4 for CDK1
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% c("CDK1", "TP53", "CASP3")) {
    if(result_list[[i]]$gene %in% c("CDK1", "TP53", "CASP3")) {
      res_low <- intersect(res_low, result_list[[i]]$low)  
    } else {
      res_low <- intersect(res_low, result_list[[i]]$high)
    }
    
    print(length(res_low))
  }
}
# [1] 237
# [1] 237
# [1] 203

# tumour classification & survival analysis
surdata$group <- 3
for (i in 1:nrow(surdata)) {
  if(row.names(surdata)[i] %in% res_high) {
    surdata$group[i] <- 1
  }
  if(row.names(surdata)[i] %in% res_low) {
    surdata$group[i] <- 2
  }
}
nrow(surdata[which(surdata$group == 1),]) # 3 genes high
nrow(surdata[which(surdata$group == 2),]) # 3 genes low
nrow(surdata[which(surdata$group == 3),]) # others
# [1] 51
# [1] 203
# [1] 94

#------------------------------------
# 3 groups
# surdata <- surdata[which(surdata$group %in% c(1,2)),]
clusterNum = length(unique(surdata$group))
time <- surdata$OS.time
status <- surdata$OS
group <- surdata$group
dataset = list(time, status, x = group)  
surv = survfit(Surv(time, status) ~ x, dataset)
mainTitle <- "Survival analysis"
if(clusterNum>1){
  sdf=NULL
  sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
  print(sdf)
  p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
}else{
  cat("There is only one cluster in the group")
  p_value=1
}

# survival chart
# draw the chart
f <- paste(outDir, "/LGG_sur_chart_3up_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
# myCol <- wes_palette("Zissou1")
myCol <- wes_palette("Darjeeling2")

# Graph 1
mar.default <- c(5,4,4,2) + 0.1 # c(bottom, left, top, right)
par(mar = mar.default + c(0, 1, 0, 0)) 
title=paste(mainTitle, " (", clusterNum, " clusters)", sep="")
plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
     main = title, font.main=2, cex.lab=1.6, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
       lty=1, lwd=3, cex=1.3, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
       seg.len = 0.3)
digit=ceiling(-log10(p_value)+2)
if (p_value < 2e-16) {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value < 2e-16"),col="blue",font=2,cex=1.5)
} else {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1.5)  
}
dev.off()

# heatmap

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
x <- surdata
x <- x[,which(colnames(x) %in% y$gene)]

# change from log count to count
x_count <- 2^x

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00
dat <- logCPM

rr <- dat
rr <- rr[which(row.names(rr) %in% c("CDK1", "TP53", "CASP3")),]
r1 <- rr[,which(colnames(rr) %in% res_high)]
ncol(r1)
r2 <- cbind(r1, rr[,which(colnames(rr) %in% res_low)])
ncol(r2) - ncol(r1)
# r3 <- cbind(r2, rr[,which(!(colnames(rr) %in% c(res_high,res_low)))])
# ncol(r3) - ncol(r2)
# r2 <- r3
# [1] 51
# [1] 203

# draw the chart
f <- paste(outDir, "/LGG_heatmap_3up_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
data <- as.matrix(r2)
my_group <- c(rep(1, 3))
rowSide <- brewer.pal(3, "Set1")[my_group]
#my_group2 <- c(rep(1, 489), rep(2,390))
my_group2 <- c(rep(1, 51), rep(2,203))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

#================================================================
# (32) Triple negative breast cancer (TNBC)
# Compare gene expression, copper genes, up genes and down genes, using count and cpm
# heatmap for clustered groups
# identify subset of genes for tumour classification & survival analysis
# run on server
#================================================================

# Get only triple negative breast cancer (i.e., Basal)

# Get cancer sub types
subtypes <- c("BRCA")
all <- PanCancerAtlas_subtypes()
PAN.path.subtypes <- all[which(all$cancer.type %in% subtypes),]
PAN.path.subtypes[1:5,1:6]
nrow(PAN.path.subtypes)
# # A tibble: 5 × 6
# pan.samplesID                cancer.type Subtype_mRNA Subtyp…¹ Subty…² Subty…³
# <chr>                        <chr>       <chr>        <chr>      <int> <chr>
#   1 TCGA-E2-A158-11A-22R-A12D-07 BRCA        Normal       NA            NA NA
# 2 TCGA-BH-A0DD-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 3 TCGA-BH-A1EO-11A-31R-A137-07 BRCA        LumA         NA            NA NA
# 4 TCGA-BH-A0B5-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 5 TCGA-A7-A13G-11A-51R-A13Q-07 BRCA        LumA         NA            NA NA
# # … with abbreviated variable names ¹​Subtype_DNAmeth, ²​Subtype_protein,
# #   ³​Subtype_miRNA
# [1] 1218
PAN.subtypes <- PAN.path.subtypes
i <- sapply(PAN.subtypes, is.factor)
PAN.subtypes[i] <- lapply(PAN.subtypes[i], as.character)
table(PAN.subtypes$Subtype_mRNA)
# Basal   Her2   LumA   LumB Normal
# 193     82    581    219    143

PAN.subtypes <- as.data.frame(PAN.subtypes)
PAN.subtypes$pan.samplesID <- substr(PAN.subtypes$pan.samplesID, 1, 15)
PAN.subtypes[1:5,1:6]
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA-E2-A158-11        BRCA       Normal            <NA>              NA
# 2 TCGA-BH-A0DD-11        BRCA         LumA            <NA>              NA
# 3 TCGA-BH-A1EO-11        BRCA         LumA            <NA>              NA
# 4 TCGA-BH-A0B5-11        BRCA         LumA            <NA>              NA
# 5 TCGA-A7-A13G-11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
PAN.subtypes$pan.samplesID <- str_replace_all(PAN.subtypes$pan.samplesID, "-", ".")
PAN.subtypes[1:5,1:6]
nrow(PAN.subtypes)
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA.E2.A158.11        BRCA       Normal            <NA>              NA
# 2 TCGA.BH.A0DD.11        BRCA         LumA            <NA>              NA
# 3 TCGA.BH.A1EO.11        BRCA         LumA            <NA>              NA
# 4 TCGA.BH.A0B5.11        BRCA         LumA            <NA>              NA
# 5 TCGA.A7.A13G.11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
#   [1] 1218

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TNBC

type <- "Basal"
patients <- PAN.subtypes[which(PAN.subtypes$Subtype_mRNA == type),1]
length(patients)
# [1] 193

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
dat <- surdata[which(rownames(surdata) %in% patients),]
nrow(dat)
dat[1:3,1:4]
# [1] 191
# cancertype OS OS.time  AANAT
# TCGA.A1.A0SK.01       BRCA  1     967 1.5850
# TCGA.A1.A0SO.01       BRCA  0     852 2.3219
# TCGA.A1.A0SP.01       BRCA  0     584 1.0000

cancertype <- "BRCA"
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y[1:4,1:5]
# gene frequency ACC  AML BLCA
# 1  CDK1        18  UP   UP   UP
# 2   ALB        13   0 DOWN    0
# 3 AP1S1        11   0 DOWN    0
# 4 CASP3         9   0    0    0
# Critical copper genes
criticalCopperGenes <- y
upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'UP'),]$gene
downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'DOWN'),]$gene
geneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] %in% c("UP", "DOWN")),]$gene
upGeneList
downGeneList
geneList
# [1] "CDK1"  "CASP3" "GSK3B" "XIAP"  "ARF1"  "SORD"
# [1] "ALB"      "MAP1LC3A" "SNCA"     "JUN"      "FOXO1"    "AOC3"     "S100A12"
# [1] "CDK1"     "ALB"      "CASP3"    "MAP1LC3A" "SNCA"     "GSK3B"
# [7] "JUN"      "XIAP"     "FOXO1"    "ARF1"     "AOC3"     "SORD"
# [13] "S100A12"

#------------------------------------
# For 3 groups
# select a set of genes which have the same pattern
# e.g., 5 genes while 3 genes high and 2 genes low in a set of samples,
# and 3 genes low and 2 genes high in another set
#------------------------------------
# For 13 genes, identify 2 groups up & down, TNBC
# Top3, including both up and down

# Get genes with high node degree
critical_nodes <- read.csv(paste(rootDir, "/Data/BRCA/critical_nodes.csv", sep = ""))
copper.critical<- critical_nodes[critical_nodes$Name %in% geneList, ]
copper.critical <- copper.critical[order(copper.critical$K, decreasing = TRUE),]
top3 <- copper.critical[1:3,]
top3
# Name   K Kin Kout TypeI TypeII
# 131  CDK1 222 180   42     0      1
# 7     ALB  73  28   45     0      1
# 68  CASP3  67  29   38     0      1
# up: CDK1, CASP3
# down: ALB

surdata <- dat

# for all 13 genes
exp_data <- surdata[,-c(1,2,3)]
exp_data <- exp_data[,colnames(exp_data) %in% geneList]
exp_data[1:4,1:5]
nrow(exp_data)
ncol(exp_data)
# ALB    AOC3    ARF1   CASP3    CDK1
# TCGA.A1.A0SK.01 11.0341  8.9769 14.1411 11.2691 13.3570
# TCGA.A1.A0SO.01  6.9773  9.7558 15.3215 10.8001 11.8665
# TCGA.A1.A0SP.01  6.0444 10.3196 14.9893 11.6874 11.4067
# TCGA.A2.A04P.01  5.8329 10.5743 15.3696  9.9469 11.5324
# [1] 191
# [1] 13

x <- exp_data

# change from log count to count
x_count <- 2^x

# compute CPM
dat <- t(x_count)
CPM <- cpm(dat, log=FALSE)

exp_data <- CPM
exp_data <- t(exp_data)
mean13 <- colMeans(exp_data)

result_list <- list()
list_name <- c(1:13)
for (i in 1:13) {
  gene <- colnames(exp_data)[i]
  high <- rownames(exp_data)[which(exp_data[,i] >= mean13[i])]
  low <- rownames(exp_data)[which(exp_data[,i] < mean13[i])]
  result_list[[list_name[i]]] <- list(gene = gene, high = high, low = low)
}

gene_set <- top3$Name
# > gene_set
# [1] "CDK1"  "ALB"   "CASP3"
# > colnames(exp_data)
# [1] "ALB"      "AOC3"     "ARF1"     "CASP3"    "CDK1"     "FOXO1"
# [7] "GSK3B"    "JUN"      "MAP1LC3A" "S100A12"  "SNCA"     "SORD"
# [13] "XIAP"
# > result_list[[4]]$gene
# [1] "CASP3"
# up: CDK1, CASP3
# down: ALB
# 2 up genes & 1 down genes
res_high <- result_list[[4]]$high # index 4 for CASP3
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% top3$Name) { # only for top 3
    if(result_list[[i]]$gene %in% c("CDK1", "CASP3")) { # 2 up genes
      res_high <- intersect(res_high, result_list[[i]]$high)  
    } else { # 1 down genes
      res_high <- intersect(res_high, result_list[[i]]$low)
    }
    print(length(res_high))
  }
}
# [1] 68
# [1] 68
# [1] 36

# 2 low & 1 high
# > result_list[[4]]$gene
# [1] "CASP3"
res_low <- result_list[[4]]$low # index 4 for CASP3
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% top3$Name) {
    if(result_list[[i]]$gene %in% c("CDK1", "CASP3")) {
      res_low <- intersect(res_low, result_list[[i]]$low)  
    } else {
      res_low <- intersect(res_low, result_list[[i]]$high)
    }
    
    print(length(res_low))
  }
}
# [1] 10
# [1] 10
# [1] 9

# tumour classification & survival analysis
surdata$group <- 3
for (i in 1:nrow(surdata)) {
  if(row.names(surdata)[i] %in% res_high) {
    surdata$group[i] <- 1
  }
  if(row.names(surdata)[i] %in% res_low) {
    surdata$group[i] <- 2
  }
}
nrow(surdata[which(surdata$group == 1),]) # 2 genes high & 1 low
nrow(surdata[which(surdata$group == 2),]) # 2 genes low & 1 high
nrow(surdata[which(surdata$group == 3),]) # others
# [1] 36
# [1] 9
# [1] 146

#------------------------------------
# 3 groups
# surdata <- surdata[which(surdata$group %in% c(1,2)),]
clusterNum = length(unique(surdata$group))
time <- surdata$OS.time
status <- surdata$OS
group <- surdata$group
dataset = list(time, status, x = group)  
surv = survfit(Surv(time, status) ~ x, dataset)
mainTitle <- "Survival analysis"
if(clusterNum>1){
  sdf=NULL
  sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
  print(sdf)
  p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
}else{
  cat("There is only one cluster in the group")
  p_value=1
}

# survival chart
# draw the chart
f <- paste(outDir, "/TNBC_sur_chart_2up1down_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
# myCol <- wes_palette("Zissou1")
myCol <- wes_palette("Darjeeling2")

# Graph 1
mar.default <- c(5,4,4,2) + 0.1 # c(bottom, left, top, right)
par(mar = mar.default + c(0, 1, 0, 0)) 
title=paste(mainTitle, " (", clusterNum, " clusters)", sep="")
plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
     main = title, font.main=2, cex.lab=1.6, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
       lty=1, lwd=3, cex=1.3, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
       seg.len = 0.3)
digit=ceiling(-log10(p_value)+2)
if (p_value < 2e-16) {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value < 2e-16"),col="blue",font=2,cex=1.5)
} else {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1.5)  
}
dev.off()

# heatmap

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
x <- surdata
x <- x[,which(colnames(x) %in% geneList)]
x <- x[which(rownames(x) %in% patients),]

# change from log count to count
x_count <- 2^x

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
# > dat[1:3,1:4]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01
# AANAT         7.000219           1.000          4.0000        7.000219
# ADAM10     3600.046934        6566.915       1924.9430      946.986468
# ADAM17      532.633736        1998.497        630.5644      301.873436
# > logCPM[1:3,1:4]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01
# AANAT         1.256488       0.3668683        1.413222       0.9366785
# ADAM10        9.469919      10.8315745        9.635874       6.9427818
# ADAM17        6.724929       9.1171051        8.029537       5.3185282
# dat <- x_count
# logCPM <- t(cpm(dat, prior.count=2, log=TRUE))
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00 # scale column (gene)
# > logCPM[1:4,1:5]
# TCGA.OR.A5J1.01 TCGA.OR.A5J2.01 TCGA.OR.A5J3.01 TCGA.OR.A5J5.01
# AANAT       -0.7507425      -1.4941859      -0.6197619      -1.0180031
# ADAM10      -1.9156365      -0.8455182      -1.7852130      -3.9017003
# ADAM17      -2.0350440      -0.5656860      -1.2337084      -2.8989047
# ALB         -0.9006157      -0.5285294      -0.4871709       0.2612733
# TCGA.OR.A5J6.01
# AANAT       -0.5880477
# ADAM10      -1.3724846
# ADAM17      -1.1933036
# ALB         -0.6570130
dat <- logCPM

rr <- dat
rr <- rr[which(row.names(rr) %in% top3$Name),]
r1 <- rr[,which(colnames(rr) %in% res_high)]
ncol(r1)
r2 <- cbind(r1, rr[,which(colnames(rr) %in% res_low)])
ncol(r2) - ncol(r1)
# r3 <- cbind(r2, rr[,which(!(colnames(rr) %in% c(res_high,res_low)))])
# ncol(r3) - ncol(r2)
# r2 <- r3
# [1] 36
# [1] 9

# draw the chart
f <- paste(outDir, "/TNBC_heatmap_2up1down_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
data <- as.matrix(r2)
my_group <- c(rep(1, 3))
rowSide <- brewer.pal(3, "Set1")[my_group]
#my_group2 <- c(rep(1, 489), rep(2,390))
my_group2 <- c(rep(1, 36), rep(2,9))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

#================================================================
# (33) Analyse the effects of cuproplasia-related genes in BRCA on survival for breast cancer subtypes
# Univariate Cox Analysis
# TNBC
#================================================================

# Get only triple negative breast cancer (i.e., Basal)

# Get cancer sub types
subtypes <- c("BRCA")
all <- PanCancerAtlas_subtypes()
PAN.path.subtypes <- all[which(all$cancer.type %in% subtypes),]
PAN.path.subtypes[1:5,1:6]
nrow(PAN.path.subtypes)
# # A tibble: 5 × 6
# pan.samplesID                cancer.type Subtype_mRNA Subtyp…¹ Subty…² Subty…³
# <chr>                        <chr>       <chr>        <chr>      <int> <chr>
#   1 TCGA-E2-A158-11A-22R-A12D-07 BRCA        Normal       NA            NA NA
# 2 TCGA-BH-A0DD-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 3 TCGA-BH-A1EO-11A-31R-A137-07 BRCA        LumA         NA            NA NA
# 4 TCGA-BH-A0B5-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 5 TCGA-A7-A13G-11A-51R-A13Q-07 BRCA        LumA         NA            NA NA
# # … with abbreviated variable names ¹​Subtype_DNAmeth, ²​Subtype_protein,
# #   ³​Subtype_miRNA
# [1] 1218
PAN.subtypes <- PAN.path.subtypes
i <- sapply(PAN.subtypes, is.factor)
PAN.subtypes[i] <- lapply(PAN.subtypes[i], as.character)
table(PAN.subtypes$Subtype_mRNA)
# Basal   Her2   LumA   LumB Normal
# 193     82    581    219    143

PAN.subtypes <- as.data.frame(PAN.subtypes)
PAN.subtypes$pan.samplesID <- substr(PAN.subtypes$pan.samplesID, 1, 15)
PAN.subtypes[1:5,1:6]
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA-E2-A158-11        BRCA       Normal            <NA>              NA
# 2 TCGA-BH-A0DD-11        BRCA         LumA            <NA>              NA
# 3 TCGA-BH-A1EO-11        BRCA         LumA            <NA>              NA
# 4 TCGA-BH-A0B5-11        BRCA         LumA            <NA>              NA
# 5 TCGA-A7-A13G-11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
PAN.subtypes$pan.samplesID <- str_replace_all(PAN.subtypes$pan.samplesID, "-", ".")
PAN.subtypes[1:5,1:6]
nrow(PAN.subtypes)
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA.E2.A158.11        BRCA       Normal            <NA>              NA
# 2 TCGA.BH.A0DD.11        BRCA         LumA            <NA>              NA
# 3 TCGA.BH.A1EO.11        BRCA         LumA            <NA>              NA
# 4 TCGA.BH.A0B5.11        BRCA         LumA            <NA>              NA
# 5 TCGA.A7.A13G.11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
#   [1] 1218

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TNBC

type <- "Basal"
patients <- PAN.subtypes[which(PAN.subtypes$Subtype_mRNA == type),1]
length(patients)
# [1] 193

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
dat <- surdata[which(rownames(surdata) %in% patients),]
nrow(dat)
dat[1:3,1:4]
# [1] 191
# cancertype OS OS.time  AANAT
# TCGA.A1.A0SK.01       BRCA  1     967 1.5850
# TCGA.A1.A0SO.01       BRCA  0     852 2.3219
# TCGA.A1.A0SP.01       BRCA  0     584 1.0000

cancertype <- "BRCA"
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y[1:4,1:5]
# gene frequency ACC  AML BLCA
# 1  CDK1        18  UP   UP   UP
# 2   ALB        13   0 DOWN    0
# 3 AP1S1        11   0 DOWN    0
# 4 CASP3         9   0    0    0
# Critical copper genes
criticalCopperGenes <- y
upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'UP'),]$gene
downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'DOWN'),]$gene
geneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] %in% c("UP", "DOWN")),]$gene
upGeneList
downGeneList
geneList
# [1] "CDK1"  "CASP3" "GSK3B" "XIAP"  "ARF1"  "SORD"
# [1] "ALB"      "MAP1LC3A" "SNCA"     "JUN"      "FOXO1"    "AOC3"     "S100A12"
# [1] "CDK1"     "ALB"      "CASP3"    "MAP1LC3A" "SNCA"     "GSK3B"
# [7] "JUN"      "XIAP"     "FOXO1"    "ARF1"     "AOC3"     "SORD"
# [13] "S100A12"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#-----------------------------------------------------------
# nrow(dat)
# dat[1:3,1:4]
# # [1] 191
# # cancertype OS OS.time  AANAT
# # TCGA.A1.A0SK.01       BRCA  1     967 1.5850
# # TCGA.A1.A0SO.01       BRCA  0     852 2.3219
# # TCGA.A1.A0SP.01       BRCA  0     584 1.0000

# Univariate COX Analysis

# to check distibution of gene expression (Optional) - they should have a normal distribution
# check gene MT1X
#ggplot(coxdata, aes(x = MT1X)) + geom_histogram(color = "black", fill = "white") 

res_mod = ezcox(dat, time = "OS.time", status = "OS", covariates = geneList, global_method = c("likelihood", "wald", "logrank"), return_models = TRUE)
mds = get_models(res_mod)
str(mds, max.level = 1)
show_models(mds)

cox_res <- ezcox(dat, time = "OS.time", status = "OS", covariates = geneList, global_method = c("likelihood", "wald", "logrank"))
unicox <- cox_res
unicox[1:5,1:6]
nrow(unicox)
ncol(unicox)
# # A tibble: 5 × 6
# Variable is_control contrast_level ref_level n_contrast n_ref
# <chr>    <lgl>      <chr>          <chr>          <int> <int>
#   1 CDK1     FALSE      CDK1           CDK1             191   191
# 2 ALB      FALSE      ALB            ALB              191   191
# 3 CASP3    FALSE      CASP3          CASP3            191   191
# 4 MAP1LC3A FALSE      MAP1LC3A       MAP1LC3A         191   191
# 5 SNCA     FALSE      SNCA           SNCA             191   191
# [1] 13
# [1] 12

write.csv(unicox,file = paste(outDir, "/", "TNBC_unicox.csv", sep =""),row.names = FALSE)

######Get significant genes#######

unicox_dif <- unicox[which(unicox$global.pval <= 0.05), ]
nrow(unicox_dif)
# > nrow(unicox_dif)
# [1] 0

# > unicox[,c(1,12)]
# # A tibble: 13 × 2
# Variable global.pval
# <chr>          <dbl>
#   1 CDK1          0.779
# 2 ALB           0.952
# 3 CASP3         0.236
# 4 MAP1LC3A      0.521
# 5 SNCA          0.97
# 6 GSK3B         0.418
# 7 JUN           0.0513
# 8 XIAP          0.503
# 9 FOXO1         0.943
# 10 ARF1          0.513
# 11 AOC3          0.363
# 12 SORD          0.825
# 13 S100A12       0.997

#-----------------------------------------------------------

#================================================================

#================================================================
# (34) BRCA
# Compare gene expression, copper genes, up genes and down genes, using count and cpm
# heatmap for clustered groups
# identify subset of genes for tumour classification & survival analysis
# run on server
#================================================================

# Get cancer sub types
subtypes <- c("BRCA")
all <- PanCancerAtlas_subtypes()
PAN.path.subtypes <- all[which(all$cancer.type %in% subtypes),]
PAN.path.subtypes[1:5,1:6]
nrow(PAN.path.subtypes)
# # A tibble: 5 × 6
# pan.samplesID                cancer.type Subtype_mRNA Subtyp…¹ Subty…² Subty…³
# <chr>                        <chr>       <chr>        <chr>      <int> <chr>
#   1 TCGA-E2-A158-11A-22R-A12D-07 BRCA        Normal       NA            NA NA
# 2 TCGA-BH-A0DD-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 3 TCGA-BH-A1EO-11A-31R-A137-07 BRCA        LumA         NA            NA NA
# 4 TCGA-BH-A0B5-11A-23R-A12P-07 BRCA        LumA         NA            NA NA
# 5 TCGA-A7-A13G-11A-51R-A13Q-07 BRCA        LumA         NA            NA NA
# # … with abbreviated variable names ¹​Subtype_DNAmeth, ²​Subtype_protein,
# #   ³​Subtype_miRNA
# [1] 1218
PAN.subtypes <- PAN.path.subtypes
i <- sapply(PAN.subtypes, is.factor)
PAN.subtypes[i] <- lapply(PAN.subtypes[i], as.character)
table(PAN.subtypes$Subtype_mRNA)
# Basal   Her2   LumA   LumB Normal
# 193     82    581    219    143

PAN.subtypes <- as.data.frame(PAN.subtypes)
PAN.subtypes$pan.samplesID <- substr(PAN.subtypes$pan.samplesID, 1, 15)
PAN.subtypes[1:5,1:6]
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA-E2-A158-11        BRCA       Normal            <NA>              NA
# 2 TCGA-BH-A0DD-11        BRCA         LumA            <NA>              NA
# 3 TCGA-BH-A1EO-11        BRCA         LumA            <NA>              NA
# 4 TCGA-BH-A0B5-11        BRCA         LumA            <NA>              NA
# 5 TCGA-A7-A13G-11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
PAN.subtypes$pan.samplesID <- str_replace_all(PAN.subtypes$pan.samplesID, "-", ".")
PAN.subtypes[1:5,1:6]
nrow(PAN.subtypes)
# pan.samplesID cancer.type Subtype_mRNA Subtype_DNAmeth Subtype_protein
# 1 TCGA.E2.A158.11        BRCA       Normal            <NA>              NA
# 2 TCGA.BH.A0DD.11        BRCA         LumA            <NA>              NA
# 3 TCGA.BH.A1EO.11        BRCA         LumA            <NA>              NA
# 4 TCGA.BH.A0B5.11        BRCA         LumA            <NA>              NA
# 5 TCGA.A7.A13G.11        BRCA         LumA            <NA>              NA
# Subtype_miRNA
# 1          <NA>
#   2          <NA>
#   3          <NA>
#   4          <NA>
#   5          <NA>
#   [1] 1218

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# BRCA

patients <- PAN.subtypes[,1]
length(patients)
# [1] 1218

outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
dat <- surdata[which(rownames(surdata) %in% patients),]
nrow(dat)
dat[1:3,1:4]
# [1] 1211
# cancertype OS OS.time AANAT
# TCGA.3C.AAAU.01       BRCA  0    4047 1.000
# TCGA.3C.AALI.01       BRCA  0    4005 3.585
# TCGA.3C.AALJ.01       BRCA  0    1474 2.000

cancertype <- "BRCA"
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
y[1:4,1:5]
# gene frequency ACC  AML BLCA
# 1  CDK1        18  UP   UP   UP
# 2   ALB        13   0 DOWN    0
# 3 AP1S1        11   0 DOWN    0
# 4 CASP3         9   0    0    0
# Critical copper genes
criticalCopperGenes <- y
upGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'UP'),]$gene
downGeneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] == 'DOWN'),]$gene
geneList <- criticalCopperGenes[which(criticalCopperGenes[,cancertype] %in% c("UP", "DOWN")),]$gene
upGeneList
downGeneList
geneList
# [1] "CDK1"  "CASP3" "GSK3B" "XIAP"  "ARF1"  "SORD"
# [1] "ALB"      "MAP1LC3A" "SNCA"     "JUN"      "FOXO1"    "AOC3"     "S100A12"
# [1] "CDK1"     "ALB"      "CASP3"    "MAP1LC3A" "SNCA"     "GSK3B"
# [7] "JUN"      "XIAP"     "FOXO1"    "ARF1"     "AOC3"     "SORD"
# [13] "S100A12"

#------------------------------------
# For 3 groups
# select a set of genes which have the same pattern
# e.g., 5 genes while 3 genes high and 2 genes low in a set of samples,
# and 3 genes low and 2 genes high in another set
#------------------------------------
# For 13 genes, identify 2 groups up & down
# Top3, including both up and down

# Get genes with high node degree
critical_nodes <- read.csv(paste(rootDir, "/Data/BRCA/critical_nodes.csv", sep = ""))
copper.critical<- critical_nodes[critical_nodes$Name %in% geneList, ]
copper.critical <- copper.critical[order(copper.critical$K, decreasing = TRUE),]
top3 <- copper.critical[1:3,]
top3
# Name   K Kin Kout TypeI TypeII
# 131  CDK1 222 180   42     0      1
# 7     ALB  73  28   45     0      1
# 68  CASP3  67  29   38     0      1
# up: CDK1, CASP3
# down: ALB

surdata <- dat

# for all 13 genes
exp_data <- surdata[,-c(1,2,3)]
exp_data <- exp_data[,colnames(exp_data) %in% geneList]
exp_data[1:4,1:5]
nrow(exp_data)
ncol(exp_data)
# ALB    AOC3    ARF1   CASP3    CDK1
# TCGA.3C.AAAU.01 5.5850  9.3187 15.0249 11.0519 11.4221
# TCGA.3C.AALI.01 5.5236  9.8588 14.4873 11.2621 10.9024
# TCGA.3C.AALJ.01 3.3219 10.3401 14.4499  9.8765 11.4752
# TCGA.3C.AALK.01 5.1293 10.1260 15.0125 10.7805 10.5536
# [1] 1211
# [1] 13

x <- exp_data

# change from log count to count
x_count <- 2^x

# compute CPM
dat <- t(x_count)
CPM <- cpm(dat, log=FALSE)

exp_data <- CPM
exp_data <- t(exp_data)
mean13 <- colMeans(exp_data)

result_list <- list()
list_name <- c(1:13)
for (i in 1:13) {
  gene <- colnames(exp_data)[i]
  high <- rownames(exp_data)[which(exp_data[,i] >= mean13[i])]
  low <- rownames(exp_data)[which(exp_data[,i] < mean13[i])]
  result_list[[list_name[i]]] <- list(gene = gene, high = high, low = low)
}

gene_set <- top3$Name
# > gene_set
# [1] "CDK1"  "ALB"   "CASP3"
# > colnames(exp_data)
# [1] "ALB"      "AOC3"     "ARF1"     "CASP3"    "CDK1"     "FOXO1"
# [7] "GSK3B"    "JUN"      "MAP1LC3A" "S100A12"  "SNCA"     "SORD"
# [13] "XIAP"
# > result_list[[4]]$gene
# [1] "CASP3"
# up: CDK1, CASP3
# down: ALB
# 2 up genes & 1 down genes
res_high <- result_list[[4]]$high # index 4 for CASP3
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% top3$Name) { # only for top 3
    if(result_list[[i]]$gene %in% c("CDK1", "CASP3")) { # 2 up genes
      res_high <- intersect(res_high, result_list[[i]]$high)  
    } else { # 1 down genes
      res_high <- intersect(res_high, result_list[[i]]$low)
    }
    print(length(res_high))
  }
}
# [1] 494
# [1] 494
# [1] 287

# 2 low & 1 high
# > result_list[[4]]$gene
# [1] "CASP3"
res_low <- result_list[[4]]$low # index 4 for CASP3
for (i in 1:length(result_list)) {
  if(result_list[[i]]$gene %in% top3$Name) {
    if(result_list[[i]]$gene %in% c("CDK1", "CASP3")) {
      res_low <- intersect(res_low, result_list[[i]]$low)  
    } else {
      res_low <- intersect(res_low, result_list[[i]]$high)
    }
    
    print(length(res_low))
  }
}
# [1] 91
# [1] 91
# [1] 85

# tumour classification & survival analysis
surdata$group <- 3
for (i in 1:nrow(surdata)) {
  if(row.names(surdata)[i] %in% res_high) {
    surdata$group[i] <- 1
  }
  if(row.names(surdata)[i] %in% res_low) {
    surdata$group[i] <- 2
  }
}
nrow(surdata[which(surdata$group == 1),]) # 2 genes high & 1 low
nrow(surdata[which(surdata$group == 2),]) # 2 genes low & 1 high
nrow(surdata[which(surdata$group == 3),]) # others
# [1] 287
# [1] 85
# [1] 839

#------------------------------------
# 3 groups
# surdata <- surdata[which(surdata$group %in% c(1,2)),]
clusterNum = length(unique(surdata$group))
time <- surdata$OS.time
status <- surdata$OS
group <- surdata$group
dataset = list(time, status, x = group)  
surv = survfit(Surv(time, status) ~ x, dataset)
mainTitle <- "Survival analysis"
if(clusterNum>1){
  sdf=NULL
  sdf=survdiff(Surv(time, status) ~ group) ##log-rank test
  cat("                                                     \n")
  cat("*****************************************************\n")
  cat(paste(mainTitle, " (Number of clusters = ", clusterNum, ")", ""))
  print(sdf)
  p_value = 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
}else{
  cat("There is only one cluster in the group")
  p_value=1
}

# survival chart
# draw the chart
f <- paste(outDir, "/BRCA_sur_chart_2up1down_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
# myCol <- wes_palette("Zissou1")
myCol <- wes_palette("Darjeeling2")

# Graph 1
mar.default <- c(5,4,4,2) + 0.1 # c(bottom, left, top, right)
par(mar = mar.default + c(0, 1, 0, 0)) 
title=paste(mainTitle, " (", clusterNum, " clusters)", sep="")
plot(surv, lty = 1, col=myCol[2:(clusterNum+1)], lwd=2, xscale=30, xlab="Survival time (Months)", ylab="Survival probability",
     main = title, font.main=2, cex.lab=1.6, cex.axis=1.5, cex.main=1.7, cex.sub=1.5)
legend(x=par("usr")[2]*0.6,y=par("usr")[4]*0.95,x.intersp=0.05,y.intersp=1.2, paste(" Subtype", 1:clusterNum),
       lty=1, lwd=3, cex=1.3, text.font=2, text.col=myCol[2:(clusterNum+1)], bty="n", col=myCol[2:(clusterNum+1)],
       seg.len = 0.3)
digit=ceiling(-log10(p_value)+2)
if (p_value < 2e-16) {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value < 2e-16"),col="blue",font=2,cex=1.5)
} else {
  text(x=par("usr")[2]*0.6,y=par("usr")[4]*0.96,paste("p-value =",round(p_value,digit)),col="blue",font=2,cex=1.5)  
}
dev.off()

# heatmap

# Using data from surdata
# gene data, log count
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/surdata.Rdata",sep="")
load(fileName) # return surdata
surdata[1:5,1:4]
nrow(surdata)
ncol(surdata)
# cancertype OS OS.time  AANAT
# TCGA.OR.A5J1.01        ACC  1    1355 2.8074
# TCGA.OR.A5J2.01        ACC  1    1677 0.0000
# TCGA.OR.A5J3.01        ACC  0    2091 2.0000
# TCGA.OR.A5J5.01        ACC  1     365 2.8074
# TCGA.OR.A5J6.01        ACC  0    2703 2.0000
# [1] 6727
# [1] 60

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)
x <- surdata
x <- x[,which(colnames(x) %in% geneList)]
x <- x[which(rownames(x) %in% patients),]

# change from log count to count
x_count <- 2^x

# compute logCPM
dat <- t(x_count)
logCPM <- cpm(dat, prior.count=2, log=TRUE)
logCPM <- t(scale(t(logCPM))) # z score, need memory, qsub -I -l select=1:ncpus=2:mem=24gb,walltime=8:00:00 # scale column (gene)
dat <- logCPM

rr <- dat
rr <- rr[which(row.names(rr) %in% top3$Name),]
r1 <- rr[,which(colnames(rr) %in% res_high)]
ncol(r1)
r2 <- cbind(r1, rr[,which(colnames(rr) %in% res_low)])
ncol(r2) - ncol(r1)
# r3 <- cbind(r2, rr[,which(!(colnames(rr) %in% c(res_high,res_low)))])
# ncol(r3) - ncol(r2)
# r2 <- r3
# [1] 287
# [1] 85

# draw the chart
f <- paste(outDir, "/BRCA_heatmap_2up1down_3groups.pdf", sep = "")
pdf(file = f, width = 6, height =  6)
data <- as.matrix(r2)
my_group <- c(rep(1, 3))
rowSide <- brewer.pal(3, "Set1")[my_group]
#my_group2 <- c(rep(1, 489), rep(2,390))
my_group2 <- c(rep(1, 287), rep(2,85))
colSide <- brewer.pal(3, "Set2")[my_group2]
# colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
# heatmap(data, Colv = NA, Rowv = NA, scale="column" , RowSideColors=colSide, col=colMain   )
heatmap(data, Colv = NA, Rowv = NA, RowSideColors=rowSide, ColSideColors=colSide , cexRow = 1.5, cexCol = 1.5, margins = c(10, 10))
dev.off()

#================================================================
# (35) Lollipop
# run on server & local
#================================================================
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Server
#-----------------------------------------
# Evaluate SNVs
# Remember to remove SNVs with VAF < 0.05

threshold <- 0.05

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get 30 genes

outDir <- paste(rootDir, "/Data/Output", sep = "")

fileName<-paste(outDir, "/criticalCopperGenesUpDown.csv",sep="")
y <- read.csv(fileName)

# 21 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList21 <- y[which(y$numUp > y$numDown),]$gene

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene

# 35 cox genes
uni_gene <- read.csv(file = paste(outDir, "/prognostic_genes.csv", sep =""))
head(uni_gene)
nrow(uni_gene)
# nrow(uni_gene)
# x
# 1     CDK1
# 2    AP1S1
# 3    CASP3
# 4 MAP1LC3A
# 5     SNCA
# 6  TMPRSS6
# [1] 35

# Cox genes
geneList18 <- geneList21[which(geneList21 %in% uni_gene$x)]
geneList12 <- geneList13[which(geneList13 %in% uni_gene$x)]

geneList30 <- c(geneList18, geneList12)
geneList30
# [1] "CDK1"     "AP1S1"    "CASP3"    "TMPRSS6"  "GSK3B"    "APP"
# [7] "COX17"    "XIAP"     "ARF1"     "GPC1"     "SORD"     "ATP7A"
# [13] "SP1"      "MT-CO1"   "SLC11A2"  "ATP6AP1"  "ADAM10"   "CP"
# [19] "MAP1LC3A" "SNCA"     "MAPT"     "JUN"      "CYP1A1"   "AOC3"
# [25] "PRNP"     "DBH"      "S100A12"  "AQP1"     "MT1X"     "XAF1"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get mutation data

subtypes <- PanCancerAtlas_subtypes()
subtypes <- unique(subtypes$cancer.type)
subtypes <- c(subtypes, "PAAD")
subtypes <- subtypes[-which(subtypes %in% c("READ", "HNSC"))]
subtypes <- subtypes[order(subtypes, decreasing = FALSE)]

SNV_types <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del",
               "In_Frame_Ins")

# geneList30

# get all mutation info
n <- length(subtypes)
for (i in 1:n) {
  cancertype <- subtypes[i]
  # Read data
  dat <- as.data.frame(fread(paste0(outDir, "/SNV/", cancertype, "_maf_data.IdTrans.deleterious.csv")))
  dat <- dat[which(dat$VAF >= threshold),]
  # Only get necessary data
  dat <- dat[,c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","HGVSp_Short")]
  if(nrow(dat) > 0) {
    dat$CancerType <- cancertype
    if(i == 1) {
      res <- dat
    } else {
      res <- rbind(res,dat)
    }  
  }
}

# 30 genes
write.csv(res[which(res$Hugo_Symbol %in% geneList30),], paste0(outDir, "/SNV/SNV30Genes_005_Lollipop.csv"), row.names=FALSE, quote=FALSE)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# Local

rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"

# load data
mutation.csv <- paste0(rootDir, "/Data/SNV/SNV30Genes_005_Lollipop.csv")

# ============================================
# read in data
#   "gene.symbol.col"    : column of gene symbol
#   "variant.class.col"  : column of variant class
#   "protein.change.col" : colum of protein change column
# ============================================
mutation.dat <- readMAF(mutation.csv,
                        gene.symbol.col = "Hugo_Symbol",
                        variant.class.col = "Variant_Classification",
                        protein.change.col = "HGVSp_Short", # "amino_acid_change",
                        sep = ",")  # column-separator of csv file

# ---------------------------
# Draw charts

geneSet <- c("ATP7A", "CP", "MAPT", "APP", "TMPRSS6")

i <- 1

curGene <- geneSet[i]

# set up chart options
plot.options <- g3Lollipop.options(
  # Chart settings
  chart.width = 600,
  chart.type = "pie",
  chart.margin = list(left = 30, right = 20, top = 20, bottom = 30),
  chart.background = "white", # "#d3d3d3",
  transition.time = 300,
  # Lollipop track settings
  lollipop.track.height = 200,
  lollipop.track.background = "white", # "#d3d3d3",
  lollipop.pop.min.size = 2,
  lollipop.pop.max.size = 8,
  lollipop.pop.info.limit = 5.5,
  lollipop.pop.info.dy = "0.24em",
  lollipop.pop.info.color = "white",
  lollipop.line.color = "#a9A9A9",
  lollipop.line.width = 1,
  lollipop.circle.color = "#ffdead",
  lollipop.circle.width = 0.4,
  lollipop.label.ratio = 2,
  lollipop.label.min.font.size = 12,
  lollipop.color.scheme = "dark2",
  highlight.text.angle = 60,
  # Domain annotation track settings
  anno.height = 24,
  anno.margin = list(top = 0, bottom = 0),
  anno.background = "white", # "#d3d3d3",
  anno.bar.fill = "#a9a9a9",
  anno.bar.margin = list(top = 4, bottom = 4),
  domain.color.scheme = "pie4",
  domain.margin = list(top = 2, bottom = 2),
  domain.text.color = "black",
  domain.text.font = "11px Serif",
  # Y-axis label
  y.axis.label = "# of mutations",
  axis.label.color = "#303030",
  axis.label.alignment = "middle",
  axis.label.font = "12px Serif",
  axis.label.dy = "-1.5em",
  y.axis.line.color = "#303030",
  y.axis.line.width = 0.5,
  y.axis.line.style = "line",
  y.max.range.ratio = 1.1,
  # Chart title settings
  title.color = "#303030",
  title.text = paste0(curGene, " gene"),
  title.font = "bold 12px Serif",
  title.alignment = "middle",
  # Chart legend settings
  legend = TRUE,
  legend.margin = list(left=20, right = 0, top = 10, bottom = 5),
  legend.interactive = TRUE,
  legend.title = "Variant classification",
  # Brush selection tool
  brush = TRUE,
  brush.selection.background = "#F8F8FF",
  brush.selection.opacity = 0.3,
  brush.border.color = "#a9a9a9",
  brush.border.width = 1,
  brush.handler.color = "#303030",
  # tooltip and zoom
  tooltip = TRUE,
  zoom = TRUE
)

g3Lollipop(mutation.dat,
           gene.symbol = curGene,
           protein.change.col = "HGVSp_Short",
           btn.style = "blue", # blue-style chart download buttons
           plot.options = plot.options,
           output.filename = curGene)
#================================================================