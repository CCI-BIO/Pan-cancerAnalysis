# Local
# Clear the environment
rm(list = ls())

# Load necessary libraries if any
library(datasets)
library(taRifx)
library(xtable)
library(TCGAbiolinks)
library(UpSetR)
library(STRINGdb)
library(writexl)
library(dplyr)
library(mygene)
library(ggplot2)
library(stringr)
library(reshape)
library(plyr)
library(scales)
library(ggpubr)
library(CancerSubtypes)
library(survival)
library(wesanderson)
library("SNFtool")
library(gridExtra)
library(openxlsx)
library(RColorBrewer)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(plotly)
library(heatmaply)
library(forcats)
library(ezcox)
library(glmnet)

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
  net <- read.table(paste(outDir, "/cancer_network_analysis.txt", sep = ""), header = FALSE, sep = "\t", dec = ".", fill = TRUE)
  t[i,1] <- cancertype
  t[i,2] <- as.numeric(net[1,3])
  t[i,3] <- as.numeric(net[5,3])
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
colnames(t) <- c('Cancer type', 'Number of CCGs', 'Critical copper genes')
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
# [1] "CDK1"     "COX17"    "DBH"      "SLC11A2"  "MAP1LC3A" "ALB"      "ANG"      "ANKRD9"   "AP1S1"    "ARF1"     "APC"      "TMPRSS6"  "GPC1"     "LOXL1"    "ATOX1"    "MAPT"    
# [17] "CASP3"    "XIAP"     "GSK3B"    "AOC3"     "S100A12"  "FOXO1"    "JUN"      "SNCA"     "SORD"     "APP"      "PRNP"     "ATP6AP1"  "F5"       "ADAM10"   "ADAM17"   "MT-CO2"  
# [33] "STEAP3"   "CYP1A1"   "CP"       "MMGT1"    "BACE1"    "AP1B1"    "TP53"     "ATP7A"    "SP1"      "CCND1"    "MT-CO1"   "PRND"     "AANAT"    "AQP1"     "MT1X"     "IL1A"    
# [49] "SPATA5"   "COMMD1"   "F8"       "SUMF1"    "XAF1"     "HEPH"     "PARK7"    "LCAT"    
# > nr # number of critical copper genes
# [1] 56
# > nc # number of cancer types
# [1] 23

x <- criticalCopperGenes
x$frequency <- rowSums(criticalCopperGenes)
x <- x[order(x$frequency, decreasing = TRUE),]
# > row.names(x)
# [1] "CDK1"     "ALB"      "AP1S1"    "CASP3"    "MAP1LC3A" "SNCA"     "TMPRSS6"  "MAPT"     "GSK3B"    "JUN"      "APP"      "CYP1A1"   "COX17"    "XIAP"     "TP53"     "FOXO1"   
# [17] "ARF1"     "GPC1"     "AOC3"     "SORD"     "PRNP"     "F5"       "ATP7A"    "SP1"      "MT-CO1"   "DBH"      "SLC11A2"  "ANG"      "S100A12"  "ATP6AP1"  "ADAM10"   "MT-CO2"  
# [33] "CP"       "BACE1"    "PRND"     "AQP1"     "MT1X"     "IL1A"     "XAF1"     "ANKRD9"   "APC"      "LOXL1"    "ATOX1"    "ADAM17"   "STEAP3"   "MMGT1"    "AP1B1"    "CCND1"   
# [49] "AANAT"    "SPATA5"   "COMMD1"   "F8"       "SUMF1"    "HEPH"     "PARK7"    "LCAT"   

print(paste(row.names(x), collapse=", "))
# > print(paste(row.names(x), collapse=", "))
# [1] "CDK1, ALB, AP1S1, CASP3, MAP1LC3A, SNCA, TMPRSS6, MAPT, GSK3B, JUN, APP, CYP1A1, COX17, XIAP, TP53, FOXO1, ARF1, GPC1, AOC3, SORD, PRNP, F5, ATP7A, SP1, MT-CO1, DBH, SLC11A2, ANG, S100A12, ATP6AP1, ADAM10, MT-CO2, CP, BACE1, PRND, AQP1, MT1X, IL1A, XAF1, ANKRD9, APC, LOXL1, ATOX1, ADAM17, STEAP3, MMGT1, AP1B1, CCND1, AANAT, SPATA5, COMMD1, F8, SUMF1, HEPH, PARK7, LCAT"

important.copper.genelist <- c('SLC31A1', 'MT2A', 'ATP7A', 'ATP7B', 'MT1X', 'ATOX1', 'COX17', 'MTH1', 'SOD1', 'GSS')
intersect(important.copper.genelist, row.names(x))
# > intersect(important.copper.genelist, row.names(x))
# [1] "ATP7A" "MT1X"  "ATOX1" "COX17"

# Only get critical copper genes in more than 1 cancer type
selectedCriticalCopperGenes <- x[which(x$frequency > 1),]
print(paste(row.names(selectedCriticalCopperGenes), collapse=", "))
nrow(selectedCriticalCopperGenes)
# > print(paste(row.names(selectedCriticalCopperGenes), collapse=", "))
# [1] "CDK1, ALB, AP1S1, CASP3, MAP1LC3A, SNCA, TMPRSS6, MAPT, GSK3B, JUN, APP, CYP1A1, COX17, XIAP, TP53, FOXO1, ARF1, GPC1, AOC3, SORD, PRNP, F5, ATP7A, SP1, MT-CO1, DBH, SLC11A2, ANG, S100A12, ATP6AP1, ADAM10, MT-CO2, CP, BACE1, PRND, AQP1, MT1X, IL1A, XAF1"
# > nrow(selectedCriticalCopperGenes)
# [1] 39

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
head(stringInteractions)
nrow(stringInteractions)
# > head(stringInteractions)
# from                   to combined_score
# 1 9606.ENSP00000000233 9606.ENSP00000324287            767
# 2 9606.ENSP00000000233 9606.ENSP00000379191            600
# 3 9606.ENSP00000000233 9606.ENSP00000355759            543
# 4 9606.ENSP00000000233 9606.ENSP00000387286            730
# 5 9606.ENSP00000000233 9606.ENSP00000331342            499
# 6 9606.ENSP00000000233 9606.ENSP00000311962            457
# > nrow(stringInteractions)
# [1] 1795134

# Convert to gene names
t <- merge(stringInteractions, string_proteins, by.x = "from", by.y = "protein_external_id")
t[1:5,1:5]
# > t[1:5,1:5]
# from                   to combined_score preferred_name protein_size
# 1 9606.ENSP00000000233 9606.ENSP00000324287            767           ARF5          180
# 2 9606.ENSP00000000233 9606.ENSP00000379191            600           ARF5          180
# 3 9606.ENSP00000000233 9606.ENSP00000355759            543           ARF5          180
# 4 9606.ENSP00000000233 9606.ENSP00000387286            730           ARF5          180
# 5 9606.ENSP00000000233 9606.ENSP00000331342            499           ARF5          180
t <- t[,c(2,3,4)]
colnames(t)[3] <- "from"
head(t)
# > head(t)
# to combined_score from
# 1 9606.ENSP00000324287            767 ARF5
# 2 9606.ENSP00000379191            600 ARF5
# 3 9606.ENSP00000355759            543 ARF5
# 4 9606.ENSP00000387286            730 ARF5
# 5 9606.ENSP00000331342            499 ARF5
# 6 9606.ENSP00000311962            457 ARF5

stringInteractions <- merge(t, string_proteins, by.x = "to", by.y = "protein_external_id")
stringInteractions[1:5,1:5]
# > stringInteractions[1:5,1:5]
# to combined_score    from preferred_name protein_size
# 1 9606.ENSP00000003084            427 CYP26B1           CFTR         1480
# 2 9606.ENSP00000003084            459    M6PR           CFTR         1480
# 3 9606.ENSP00000003084            459    M6PR           CFTR         1480
# 4 9606.ENSP00000003084            427 CYP26B1           CFTR         1480
# 5 9606.ENSP00000005226            690    CFTR          USH1C          899
stringInteractions <- stringInteractions[,c(2,3,4)]
colnames(stringInteractions)[3] <- "to"
head(stringInteractions)
# > head(stringInteractions)
# combined_score    from    to
# 1            427 CYP26B1  CFTR
# 2            459    M6PR  CFTR
# 3            459    M6PR  CFTR
# 4            427 CYP26B1  CFTR
# 5            690    CFTR USH1C
# 6            690    CFTR USH1C
stringInteractions <- stringInteractions %>% dplyr::select(to, everything())
stringInteractions <- stringInteractions %>% dplyr::select(from, everything())
head(stringInteractions)
# > head(stringInteractions)
# from    to combined_score
# 1 CYP26B1  CFTR            427
# 2    M6PR  CFTR            459
# 3    M6PR  CFTR            459
# 4 CYP26B1  CFTR            427
# 5    CFTR USH1C            690
# 6    CFTR USH1C            690

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
# Enrichment analysis
outDir <- paste(rootDir, "/Data/Output", sep = "")
fileName<-paste(outDir, "/criticalCoperGenes.csv",sep="")
selectedCriticalCopperGenes <- read.csv(fileName)
cat(selectedCriticalCopperGenes$gene, sep = '\n')

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
write.csv(GO_process, file = paste(outDir, "/Enrich/KEGG_2021_Human_table.csv",
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
# [1] "Ferroptosis"                                                "Alzheimer disease"                                         
# [3] "Colorectal cancer"                                          "Mineral absorption"                                        
# [5] "Epithelial cell signaling in Helicobacter pylori infection" "p53 signaling pathway"                                     
# [7] "Measles"                                                    "Breast cancer"                                             
# [9] "Pathways of neurodegeneration"                              "AGE-RAGE signaling pathway in diabetic complications"      
# [11] "Pathways in cancer"                                         "Endometrial cancer"                                        
# [13] "Lysosome"                                                   "Human papillomavirus infection"                            
# [15] "Cellular senescence"                                        "Wnt signaling pathway"                                     
# [17] "Small cell lung cancer"                                     "Prostate cancer"                                           
# [19] "Kaposi sarcoma-associated herpesvirus infection"            "Viral carcinogenesis"         
# for (i in 1:top) {
#   if(nchar(dat$Term[i]) > 23) {
#     dat$Term[i] <- paste(substr(dat$Term[i],1,20), "...", sep = "")
#   }
# }
dat$Term[5] <- "Epithelial cell signaling"
dat$Term[10] <- "AGE-RAGE signaling pathway"
dat$Term[19] <- "Herpesvirus infection"
dat$Term
# > dat$Term
# [1] "Ferroptosis"                    "Alzheimer disease"              "Colorectal cancer"              "Mineral absorption"            
# [5] "Epithelial cell signaling"      "p53 signaling pathway"          "Measles"                        "Breast cancer"                 
# [9] "Pathways of neurodegeneration"  "AGE-RAGE signaling pathway"     "Pathways in cancer"             "Endometrial cancer"            
# [13] "Lysosome"                       "Human papillomavirus infection" "Cellular senescence"            "Wnt signaling pathway"         
# [17] "Small cell lung cancer"         "Prostate cancer"                "Herpesvirus infection"          "Viral carcinogenesis"
l <- dat$Term
drawClustergram2(dat, "", l)
dev.off()
# Table
GO_process[,4] <- gsub("[;]","; ",GO_process[,4])
print(latex.table.by(GO_process[1:20,], digits = c(0,0,0,-3,0)), include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = force)

# Reactome
top <- 20
GO_process <- read.table(paste(outDir, "/Enrich/Reactome_2016_table.txt", sep=""),
                         as.is = TRUE, sep = "\t", header = TRUE, quote="")
GO_process <- GO_process[, c(1,2,4,9)]
colnames(GO_process) <- c("Term", "Overlap", "Adjusted p-value", "Genes")
GO_process <- GO_process[which(GO_process[,3] < 0.05),]
GO_process <- GO_process[order(GO_process[,3], decreasing = FALSE),]
write.csv(GO_process, file = paste(outDir, "/Enrich/Reactome_2016_table.csv",
                                   sep=""), row.names = FALSE, quote=TRUE)
r <- prepareDataForClustergram(GO_process, termTop = top, geneTop  = 0)
rr <- t(r)
rr <- cbind(rownames(rr),rr)
colnames(rr) <- rr[1,]
rr <- rr[-1,]
colnames(rr)[1] <- "Term"
rr[1:top,1] <- GO_process$Term[1:top]
write.csv(rr, file = paste(outDir, "/Enrich/Reactome_2016_table.csv",
                           sep=""), row.names = FALSE, quote=TRUE)

# Output
# Image
dat <- read.csv(file = paste(outDir, "/Enrich/Reactome_2016_table.csv", sep=""))
f <- paste(outDir, "/Enrich/Reactome_2016_table.pdf", sep = "")
pdf(file = f,  width = 12, height = 8, onefile=FALSE)
dat$Term
# > dat$Term
# [1] "Iron uptake and transport Homo sapiens R-HSA-917937"                                                                     
# [2] "Lysosome Vesicle Biogenesis Homo sapiens R-HSA-432720"                                                                   
# [3] "Cellular responses to stress Homo sapiens R-HSA-2262752"                                                                 
# [4] "Signaling by NOTCH Homo sapiens R-HSA-157118"                                                                            
# [5] "Nef-mediates down modulation of cell surface receptors by recruiting them to clathrin adapters Homo sapiens R-HSA-164938"
# [6] "Disease Homo sapiens R-HSA-1643685"                                                                                      
# [7] "Transmembrane transport of small molecules Homo sapiens R-HSA-382551"                                                    
# [8] "trans-Golgi Network Vesicle Budding Homo sapiens R-HSA-199992"                                                           
# [9] "Clathrin derived vesicle budding Homo sapiens R-HSA-421837"                                                              
# [10] "Metal ion SLC transporters Homo sapiens R-HSA-425410"                                                                    
# [11] "The role of Nef in HIV-1 replication and disease pathogenesis Homo sapiens R-HSA-164952"                                 
# [12] "SMAC binds to IAPs Homo sapiens R-HSA-111463"                                                                            
# [13] "SMAC-mediated dissociation of IAP:caspase complexes Homo sapiens R-HSA-111464"                                           
# [14] "SMAC-mediated apoptotic response Homo sapiens R-HSA-111469"                                                              
# [15] "Apoptosis Homo sapiens R-HSA-109581"                                                                                     
# [16] "Programmed Cell Death Homo sapiens R-HSA-5357801"                                                                        
# [17] "Apoptotic cleavage of cellular proteins Homo sapiens R-HSA-111465"                                                       
# [18] "Apoptotic factor-mediated response Homo sapiens R-HSA-111471"                                                            
# [19] "Signaling by NOTCH1 t(7;9)(NOTCH1:M1580 K2555) Translocation Mutant Homo sapiens R-HSA-2660825"                          
# [20] "Constitutive Signaling by NOTCH1 t(7;9)(NOTCH1:M1580 K2555) Translocation Mutant Homo sapiens R-HSA-2660826"   
for (i in 1:top) {
  strL <- nchar(dat$Term[i])
  dat$Term[i] <- str_trim(substr(dat$Term[i],strL-12,strL))
}
dat$Term
# > dat$Term
# [1] "R-HSA-917937"  "R-HSA-432720"  "R-HSA-2262752" "R-HSA-157118"  "R-HSA-164938"  "R-HSA-1643685" "R-HSA-382551"  "R-HSA-199992"  "R-HSA-421837" 
# [10] "R-HSA-425410"  "R-HSA-164952"  "R-HSA-111463"  "R-HSA-111464"  "R-HSA-111469"  "R-HSA-109581"  "R-HSA-5357801" "R-HSA-111465"  "R-HSA-111471" 
# [19] "R-HSA-2660825" "R-HSA-2660826"
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
# [1] 56

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
# [1] "Up critical copper genes: 31"
# > print(paste(r[which(r$numUp > r$numDown),"gene"], collapse=", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, TP53, ARF1, GPC1, SORD, F5, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP, APC, ATOX1, ADAM17, STEAP3, MMGT1, AP1B1, CCND1, SPATA5, COMMD1, SUMF1, PARK7"
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
# > r <- r[which(r$frequency > 1),]
# > print(paste("Up critical copper genes: ", nrow(r[which(r$numUp > r$numDown),]), sep = ""))
# [1] "Up critical copper genes: 20"
# > print(paste(r[which(r$numUp > r$numDown),"gene"], collapse=", "))
# [1] "CDK1, AP1S1, CASP3, TMPRSS6, GSK3B, APP, COX17, XIAP, TP53, ARF1, GPC1, SORD, F5, ATP7A, SP1, MT-CO1, SLC11A2, ATP6AP1, ADAM10, CP"
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
surdata <- obtainSurData(geneList) # save in surdata.Rdata for 56 critical copper genes
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
# [1] 56

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
dat <- dat[1:10,]
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
pdf(file = outFile,  width = 9, height = 3, onefile=FALSE)
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
# [1] 59

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
# [1] 56
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
# [1] 56
# > ncol(y)
# [1] 27
# > colnames(y)
# [1] "gene"      "frequency" "ACC"       "AML"       "BLCA"      "BRCA"
# [7] "COAD"      "ESCA"      "GBM"       "KICH"      "KIRC"      "KIRP"
# [13] "LGG"       "LIHC"      "LUAD"      "LUSC"      "OVCA"      "PAAD"
# [19] "PCPG"      "PRAD"      "SKCM"      "STAD"      "THCA"      "UCEC"
# [25] "UCS"       "numUp"     "numDown"

# All 56 critical copper genes
geneList56 <- x$gene
f56 <- paste(outDir, "/uni56.csv", sep = "")

# 31 up critical copper genes
geneList31 <- y[which(y$numUp > y$numDown),]$gene
f31 <- paste(outDir, "/uni31.csv", sep = "")

# 19 down critical copper genes
geneList19 <- y[which(y$numUp < y$numDown),]$gene
f19 <- paste(outDir, "/uni19.csv", sep = "")

# 39 popular critical cooper genes
geneList39 <- x[which(x$frequency > 1),]$gene
f39 <- paste(outDir, "/uni39.csv", sep = "")

# 20 up critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList20 <- y[which(y$numUp > y$numDown),]$gene
f20 <- paste(outDir, "/uni20.csv", sep = "")

# 13 down critical copper genes (more than 1 cancer type)
y <- y[which(y$frequency > 1),]
geneList13 <- y[which(y$numUp < y$numDown),]$gene
f13 <- paste(outDir, "/uni13.csv", sep = "")

# Univariate Cox Analysis must use log2(count+1) data

# A. First select copper genes from critical list

# geneList39
# geneList20
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
# [1] 59

# to check distibution of gene expression (Optional) - they should have a normal distribution
# check gene MT1X
#ggplot(coxdata, aes(x = MT1X)) + geom_histogram(color = "black", fill = "white") 

res_mod = ezcox(coxdata, time = "OS.time", status = "OS", covariates = geneList39, global_method = c("likelihood", "wald", "logrank"), return_models = TRUE) # for TCGA test data
mds = get_models(res_mod)
str(mds, max.level = 1)
show_models(mds)

cox_res <- ezcox(coxdata, time = "OS.time", status = "OS", covariates = geneList39, global_method = c("likelihood", "wald", "logrank"))
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
# [1] 39
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
# [1] 39
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
