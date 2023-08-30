# Katana
# Clear the environment
rm(list = ls())

# Load necessary libraries if any
# Katana
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
# AML
# ESCA
# KICH
# COAD
# STAD
# KIRC
# KIRP
# LIHC
# LUAD
# LUSC
# OVCA
# PCPG
# PRAD
# SKCM
# THCA
# UCEC
# UCS
# READ
# HNSC
# PAAD
# dir.create(file.path("C:/Users/vpham/Documents/002NetworkAnalysis/Data", "BRCA"), showWarnings = FALSE)
# rootDir="C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# outDir="C:/Users/vpham/Documents/002NetworkAnalysis/Data/LGG"
# controlDir="C:/Users/vpham/Documents/002NetworkAnalysis/control"
# cancertype="LGG"
# primarySite=""
# # Gadi
# dir.create(file.path("/scratch/eu82/vp8928/002NetworkAnalysis/Data", "AML"), showWarnings = FALSE)
# rootDir="/scratch/eu82/vp8928/002NetworkAnalysis" # And put the input files in "rootDir/Data"
# outDir="/scratch/eu82/vp8928/002NetworkAnalysis/Data/AML"
# controlDir="/scratch/eu82/vp8928/002NetworkAnalysis/control"
# cancertype="AML"
# primarySite="Blood"
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
# (3) Identifying critical nodes (Katana)
#================================================================

# Analyse controllability of the network
# Read the network
interactions <- read.csv(paste(outDir, "/cancer_network.csv",
                               sep = ""))
# Write the edges of the network for analysing controllability
dir.create(file.path(outDir, "Controllability"), showWarnings = FALSE)
write.table(interactions, paste(outDir, "/Controllability/edges.dat", sep = ""),
            row.names = FALSE, col.names=FALSE, quote=FALSE)
# Run the controllability analysis
# cmd <- paste(controlDir, "/Parse.exe ", outDir, 
#              "/Controllability/edges.dat", sep = "")
# For Linux
cmd <- paste(controlDir, "/Parse ", outDir,
             "/Controllability/edges.dat", sep = "")
system(cmd)
# cmd <- paste(controlDir, "/ControllabilityAnalysis.exe ", outDir, 
#              "/Controllability/edges.dat", sep = "")
# For Linux
cmd <- paste(controlDir, "/ControllabilityAnalysis ", outDir,
             "/Controllability/edges.dat", sep = "")
system(cmd)
# Analyse controllability of the network and output in a file
analyseControllability(paste(outDir, "/Controllability/edges.dat.output", sep = ""),
                       paste(outDir, "/analyseControllability.txt", sep = ""))

# Identify critical nodes in the network
# Read the result
nodetype <- read.table(paste(outDir, "/Controllability/edges.dat.nodetype", sep = ""))
colnames(nodetype) <- c("Name", "K", "Kin", "Kout", "TypeI", "TypeII")
# Critical nodes of the network
critical_nodes <- nodetype[which(nodetype$TypeI == 0),]
# Save file
write.csv(critical_nodes, paste(outDir, "/critical_nodes.csv", sep = ""),
          row.names = FALSE)

f <- paste(outDir, "/", cancertype, "_output.txt", sep = "")
write("Number of critical nodes",file=f,append=TRUE)
write(paste(nrow(critical_nodes), " nodes", sep = ""),file=f,append=TRUE)

#================================================================
# (4) Univariate Cox Analysis (Katana)
#================================================================

# Univariate Cox Analysis must use log2(count+1) data

# A. First select copper genes from critical list
# copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'CA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')
copper.genelist <- c('ABCB6', 'ANKRD9', 'SLC31A1', 'SLC31A2', 'PRND', 'CCDC22', 'APP', 'ARF1', 'MT2A', 'ATOX1', 'ATP7A', 'ATP7B', 'PRNP', 'SCO1', 'COX19', 'SCO2', 'CYP1A1', 'DAXX', 'BACE1', 'AOC1', 'MT1DP', 'HSF1', 'AQP1', 'AQP2', 'MT1A', 'MT1B', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1M', 'MT1X', 'MT3', 'NFE2L2', 'MT1HL1', 'SNCA', 'MAP1LC3A', 'MT4', 'BECN1', 'COMMD1', 'XIAP', 'CUTC', 'STEAP2', 'STEAP3', 'STEAP4', 'SLC11A2', 'COX17', 'CP', 'FKBP4', 'HEPHL1', 'MMGT1', 'HEPH', 'PARK7', 'AANAT', 'IL1A', 'LCAT', 'LOXL2', 'MT-CO1', 'PAM', 'ATP5F1D', 'SOD1', 'SOD3', 'SORD', 'TFRC', 'CDK1', 'MOXD2P', 'MTCO2P12', 'COX11', 'LACC1', 'DBH', 'DCT', 'ALB', 'F5', 'F8', 'OR5AR1', 'ADNP', 'ATP13A2', 'MOXD1', 'GPC1', 'ANG', 'SUMF1', 'AOC2', 'SNAI3', 'APOA4', 'COA6', 'LOX', 'LOXL1', 'MT-CO2', 'ACR', 'P2RX4', 'CUTA', 'HAMP', 'S100A5', 'S100A12', 'S100A13', 'SNCB', 'SNCG', 'TP53', 'TYR', 'LOXL4', 'LOXL3', 'AOC3', 'RNF7', 'CCS', 'AP1S1', 'AP1B1', 'TMPRSS6', 'SPATA5', 'COG2', 'ATP6V0A2', 'ATP6AP1', 'ADAM10', 'AKT1', 'MTF2', 'FOXO1', 'FOXO3', 'STEAP1', 'GSK3B', 'APC', 'JUN', 'MAPT', 'MDM2', 'MT1JP', 'MT1L', 'MTF1', 'PIK3CA', 'XAF1', 'PTEN', 'CCND1', 'SP1', 'ADAM17', 'CASP3', 'ADAM9')
length(copper.genelist)
# > length(copper.genelist)
# [1] 133

# important.copper.genelist <- c('SLC31A1', 'MT2A', 'ATP7A', 'ATP7B', 'MT1X', 'ATOX1', 'COX17', 'MTH1', 'SOD1', 'GSS')

# getting the subset DataFrame after checking values if belonging to vector
critical_nodes <- read.csv(paste(outDir, "/critical_nodes.csv", sep = ""))
copper.critical<- critical_nodes[critical_nodes$Name %in% copper.genelist, ]
copper.critical <- c(copper.critical$Name)
copper.critical
length(copper.critical)
# > copper.critical
# [1] "MAPT"   "MDM2"   "APP"    "ADAM10" "SNCA"   "GSK3B"  "AKT1"   "CASP3" 
# [9] "PTEN"   "FOXO1"  "JUN"    "SP1"    "TP53"   "DAXX"   "CP"    
# > length(copper.critical)
# [1] 15

f <- paste(outDir, "/", cancertype, "_output.txt", sep = "")
write("In copper gene list",file=f,append=TRUE)
write(paste(length(copper.critical), " genes", sep = ""),file=f,append=TRUE)
write(copper.critical,file=f,append=TRUE)

# B. Prepare PAN.log2.counts dataframe

# save(PAN.log2.counts.unique, file=paste(outDir, '/PAN.log2.counts.unique.Rdata', sep = ""))
load(paste(outDir, '/PAN.log2.counts.unique.Rdata', sep = ""))  

# create a subset of PAN.log2.counts that contains only the copper.critical
PAN.log2.copper.critical <- PAN.log2.counts.unique[PAN.log2.counts.unique$gene %in% copper.critical, ]
PAN.log2.copper.critical[1:5,1:3]
nrow(PAN.log2.copper.critical)
ncol(PAN.log2.copper.critical)
# > PAN.log2.copper.critical[1:5,1:3]
# gene TCGA.19.1787.01 TCGA.S9.A7J2.01
# 4429  ACR          2.5850          4.0875
# 5893  ALB          3.1699          5.9542
# 6234 AOC1          2.0000          1.0000
# 6235 AOC2          4.8074          5.2095
# 6236 AOC3          7.1690          6.6294
# > nrow(PAN.log2.copper.critical)
# [1] 31
# > ncol(PAN.log2.copper.critical)
# [1] 10531

# C. Need to recreate the coxdata table

if(cancertype != "PAAD") {
  x <- paste(rootDir, "/Data/Survival_SupplementalTable_S1_20171025_xena_sp", sep ="")
  PAN.surv <- read.table(x, header = TRUE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t")
  PAN.surv <- PAN.surv[, c("sample","OS", "X_PATIENT", "OS.time")]
  head(PAN.surv)
  nrow(PAN.surv)
  # > head(PAN.surv)
  # sample OS    X_PATIENT OS.time
  # 1 TCGA-OR-A5J1-01  1 TCGA-OR-A5J1    1355
  # 2 TCGA-OR-A5J2-01  1 TCGA-OR-A5J2    1677
  # 3 TCGA-OR-A5J3-01  0 TCGA-OR-A5J3    2091
  # 4 TCGA-OR-A5J4-01  1 TCGA-OR-A5J4     423
  # 5 TCGA-OR-A5J5-01  1 TCGA-OR-A5J5     365
  # 6 TCGA-OR-A5J6-01  0 TCGA-OR-A5J6    2703
  # > nrow(PAN.surv)
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
#PAN.log2.copper.critical.T[1:5,1:6]
nrow(PAN.log2.copper.critical.T)
ncol(PAN.log2.copper.critical.T)
# > PAN.log2.copper.critical.T[1:5,1:6]
# ACR     ALB    AOC1    AOC2    AOC3   AP1S1
# TCGA.19.1787.01  2.5850  3.1699  2.0000  4.8074  7.1690 10.4979
# TCGA.S9.A7J2.01  4.0875  5.9542  1.0000  5.2095  6.6294 11.6511
# TCGA.G3.A3CH.11  4.0000 22.8258  2.5850  4.3923  9.8443 10.0980
# TCGA.EK.A2RE.01  2.0000  0.0000  8.4429  5.0444  7.3219 10.5157
# TCGA.44.6778.01  3.2373  2.0000  6.4094  6.0661 12.7176 10.4878
# > nrow(PAN.log2.copper.critical.T)
# [1] 10530
# > ncol(PAN.log2.copper.critical.T)
# [1] 31

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
# > head(PAN.surv)
# sample OS    X_PATIENT OS.time       sample_ID
# 1 TCGA-OR-A5J1-01  1 TCGA-OR-A5J1    1355 TCGA.OR.A5J1.01
# 2 TCGA-OR-A5J2-01  1 TCGA-OR-A5J2    1677 TCGA.OR.A5J2.01
# 3 TCGA-OR-A5J3-01  0 TCGA-OR-A5J3    2091 TCGA.OR.A5J3.01
# 4 TCGA-OR-A5J4-01  1 TCGA-OR-A5J4     423 TCGA.OR.A5J4.01
# 5 TCGA-OR-A5J5-01  1 TCGA-OR-A5J5     365 TCGA.OR.A5J5.01
# 6 TCGA-OR-A5J6-01  0 TCGA-OR-A5J6    2703 TCGA.OR.A5J6.01
PAN.surv.tmp <- PAN.surv[ , -c(1,3)]
head(PAN.surv.tmp)
# > head(PAN.surv.tmp)
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
# PAN.log2.coxdata[1:5,1:6]
# nrow(PAN.log2.coxdata)
# ncol(PAN.log2.coxdata)
# > PAN.log2.coxdata[1:5,1:6]
# sample OS OS.time    ACR    ALB   AOC1
# 1 TCGA.02.0047.01  1     448 3.5850 3.5850 1.0000
# 2 TCGA.02.0055.01  1      76 2.1731 5.5850 2.8074
# 3 TCGA.02.2483.01  0     466 4.4594 2.5850 1.5850
# 4 TCGA.02.2485.01  0     470 2.5850 4.3219 2.3219
# 5 TCGA.04.1331.01  1    1336 5.2180 5.6147 8.6511
# > nrow(PAN.log2.coxdata)
# [1] 10491
# > ncol(PAN.log2.coxdata)
# [1] 34

#convert first column to rowname
rownames(PAN.log2.coxdata) <- PAN.log2.coxdata[,1]
PAN.log2.coxdata[,1] <- NULL
# PAN.log2.coxdata[1:5,1:6]
# > PAN.log2.coxdata[1:5,1:6]
# OS OS.time    ACR    ALB   AOC1   AOC2
# TCGA.02.0047.01  1     448 3.5850 3.5850 1.0000 3.7004
# TCGA.02.0055.01  1      76 2.1731 5.5850 2.8074 5.8074
# TCGA.02.2483.01  0     466 4.4594 2.5850 1.5850 5.6147
# TCGA.02.2485.01  0     470 2.5850 4.3219 2.3219 5.9307
# TCGA.04.1331.01  1    1336 5.2180 5.6147 8.6511 5.7004

####----First Univariate COX analysis to identify genes to be included in LASSO Regression Model----####

# Univariate COX Analysis

coxdata <- as.data.frame(PAN.log2.coxdata) # TCGA test data

# to check distibution of gene expression (Optional) - they should have a normal distribution
# check gene MT1X
#ggplot(coxdata, aes(x = MT1X)) + geom_histogram(color = "black", fill = "white") 

res_mod = ezcox(coxdata, time = "OS.time", status = "OS", covariates = copper.critical, global_method = c("likelihood", "wald", "logrank"), return_models = TRUE) # for TCGA test data
mds = get_models(res_mod)
str(mds, max.level = 1)
show_models(mds)

cox_res <- ezcox(coxdata, time = "OS.time", status = "OS", covariates = copper.critical, global_method = c("likelihood", "wald", "logrank"))
unicox <- cox_res
unicox[1:5,1:6]
nrow(unicox)
ncol(unicox)
# > unicox[1:5,1:6]
# # A tibble: 5 x 6
# Variable is_control contrast_level ref_level n_contrast n_ref
# <chr>    <lgl>      <chr>          <chr>          <int> <int>
#   1 ACR      FALSE      ACR            ACR            10429 10429
# 2 ALB      FALSE      ALB            ALB            10429 10429
# 3 AOC1     FALSE      AOC1           AOC1           10429 10429
# 4 AOC2     FALSE      AOC2           AOC2           10429 10429
# 5 AOC3     FALSE      AOC3           AOC3           10429 10429
# > nrow(unicox)
# [1] 31
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
#   Variable is_control contrast_level ref_level n_contrast n_ref
#   <chr>    <lgl>      <chr>          <chr>          <int> <int>
# 1 ACR      FALSE      ACR            ACR            10429 10429
# 2 ALB      FALSE      ALB            ALB            10429 10429
# 3 AOC1     FALSE      AOC1           AOC1           10429 10429
# 4 AOC2     FALSE      AOC2           AOC2           10429 10429
# 5 AOC3     FALSE      AOC3           AOC3           10429 10429
# > nrow(unicox)
# [1] 31
# > ncol(unicox)
# [1] 12

unicox_dif <- unicox[which(unicox$global.pval < 0.05), ]
nrow(unicox_dif)
# > nrow(unicox_dif)
# [1] 23
uni_gene <- as.character(unicox_dif$Variable)
uni_gene
# > uni_gene
#  [1] "ACR"     "ALB"     "AOC3"    "AP1S1"   "APP"     "AQP2"    "CA6"     "CCND1"   "CDK1"    "CYP1A1"  "DBH"     "DCT"     "HAMP"    "MAPT"   
# [15] "MT1A"    "MT1JP"   "MT1L"    "PRND"    "S100A12" "SNAI3"   "SOD3"    "STEAP4"  "TMPRSS6"

f <- paste(outDir, "/", cancertype, "_output.txt", sep = "")
write("Univariate Cox analysis",file=f,append=TRUE)
write(paste(nrow(unicox_dif), " genes", sep = ""),file=f,append=TRUE)
write(uni_gene,file=f,append=TRUE)

write.csv(uni_gene,file = paste(outDir, "/", cancertype, "_prognostic_genes.csv", sep =""),row.names = FALSE)

