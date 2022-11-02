# Clear the environment
rm(list = ls())

# Libraries
library(edgeR)
library(gtools)
library(gplots)
library(heatmap.plus)
library(SPIA)
library(ggfortify)
library(clusterProfiler)
library(fgsea)
library(mygene)
library(tidyverse)
library(RDAVIDWebService)
library(stats)
library(data.table)

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

#---------------------------------------
# Set environment variables if any
# Please remember to create necessary folders
# Local
rootDir <- "C:/Users/vpham/Documents/002NetworkAnalysis" # And put the input files in "rootDir/data"
outDir <- paste(rootDir, "/Data/BRCA/RNASeq", sep = "")
# geneList <- read.table(paste(rootDir, "/data/geneListErin.txt", sep = ""), header = FALSE, sep = "")
geneList <- c("ALB", "ARF1", "CASP3", "XIAP", "GSK3B", "AOC3", "S100A12", "CDK1", "FOXO1", "JUN", "MAP1LC3A", "SNCA", "SORD")
registeredDAVIDemail <- "vpham@ccia.org.au" # Change email
FC_threshold <- 2
# FC_threshold <- 1.5
GMTdirLoc <- "R:/KCA/Resources/MSigDB/" # Directory containing all MSigDB .gmt files

pathway1 <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
# pathway2 <- "VERRECCHIA_EARLY_RESPONSE_TO_TGFB1"
pathway2 <- "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN"
pathway3 <- "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP"
checkedPathways <- c(pathway1, pathway2, pathway3)
#---------------------------------------

# Set working directory
setwd(paste(rootDir, "/Data/BRCA", sep = ""))

# Include the script of functions
# source(paste(rootDir, "/script/util.R", sep=""))

#================================================================
# Analysis
#================================================================

# # VIRMA (KIAA1429)
# # HAKAI (CBLL1) duplicated
# # HNRNPG (RBMX) 
# # NSUN1 (NOP2) 
# # DNMT2 (TRDMT1) duplicated
# geneList <- rbind(geneList, "KIAA1429",  "RBMX", "NOP2")

# tpm<-read.delim("GeneExpression_TPM_Counts.txt",sep="\t",header=T)
#counts<-read.delim(paste(rootDir, "/Data/BRCA/featureCounts/featureCounts.txt", sep = ""),sep="\t",header=T)
counts <- read.table(paste(rootDir, "/Data/BRCA/featureCounts/featureCounts.txt", sep = ""), skip = 1, header = TRUE, sep ='\t')

# # Remove transcript.IDs columns if they exist
# tpm <- tpm[,-2]
# counts <- counts[,-2]

print(paste("Number of checked genes: ", length(geneList), sep = ""))
print(paste("Number of checked genes in the dataset: ", nrow(counts[which(counts$Geneid %in% geneList),]), sep = ""))
if(nrow(counts[which(counts$Geneid %in% geneList),]) < length(geneList)) {
  print("Genes are not in the dataset:")
  print(geneList[which(!(geneList %in% counts$Geneid))])
} else {
  print("All genes are in the dataset")
}
# > print(paste("Number of checked genes: ", length(geneList), sep = ""))
# [1] "Number of checked genes: 13"
# > print(paste("Number of checked genes in the dataset: ", nrow(counts[which(counts$Geneid %in% geneList),]), sep = ""))
# [1] "Number of checked genes in the dataset: 13"
# > if(nrow(counts[which(counts$Geneid %in% geneList),]) < length(geneList)) {
#   +     print("Genes are not in the dataset:")
#   +     print(geneList[which(!(geneList %in% counts$Geneid))])
#   + } else {
#     +     print("All genes are in the dataset")
#     + }
# [1] "All genes are in the dataset" 

write.table(counts[which(counts$Geneid %in% geneList),], paste(outDir, "/GeneExpression_rawCounts_list.txt", sep=""), append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
# write.table(tpm[which(tpm$gene_id %in% geneList),], paste(outDir, "/GeneExpression_TPM_Counts_list.txt", sep=""), append = FALSE, sep = " ", dec = ".",
#             row.names = TRUE, col.names = TRUE)

geneList <- geneList[which(geneList %in% counts$Geneid)]

counts <- counts[,c(1,19:ncol(counts))]
head(counts)
# > head(counts)
# Geneid MB231_CON_1_S17_L001.sorted.bam MB231_CON_2_S7_L001.sorted.bam MB231_CON_3_S13_L001.sorted.bam MB231_TEPA24_1_S8_L001.sorted.bam
# 1     DDX11L1                               0                              0                               0                                 1
# 2      WASH7P                              76                             97                              72                                97
# 3   MIR6859-1                               5                             11                               6                                13
# 4 MIR1302-2HG                               3                              0                               0                                 0
# 5   MIR1302-2                               0                              0                               0                                 0
# 6     FAM138A                               0                              0                               0                                 0
# MB231_TEPA24_2_S4_L001.sorted.bam MB231_TEPA24_3_S11_L001.sorted.bam MB231_TEPA8_1_S2_L001.sorted.bam MB231_TEPA8_2_S10_L001.sorted.bam
# 1                                 0                                  0                                0                                 0
# 2                                80                                 76                              192                                90
# 3                                12                                 16                               24                                17
# 4                                 0                                  0                                0                                 0
# 5                                 0                                  0                                0                                 0
# 6                                 0                                  0                                0                                 0
# MB231_TEPA8_3_S14_L001.sorted.bam
# 1                                 0
# 2                               153
# 3                                22
# 4                                 0
# 5                                 0
# 6                                 0                             

keep<-rowSums(cpm(counts[,2:ncol(counts)])>1)>=(0.5*ncol(counts[,2:ncol(counts)])+1)
counts<-counts[keep,]
# tpm<-tpm[keep,]
# colnames(tpm)[1]<-"genes"

print(paste("Number of checked genes: ", length(geneList), sep = ""))
print(paste("Number of checked genes after removing genes with low expression: ", nrow(counts[which(counts$Geneid %in% geneList),]), sep = ""))
if(nrow(counts[which(counts$Geneid %in% geneList),]) < length(geneList)) {
  print("Genes are removed due to low expression:")
  print(geneList[which(!(geneList %in% counts$Geneid))])
} else {
  print("All genes are kept after removing genes with low expression")
}
# > print(paste("Number of checked genes: ", length(geneList), sep = ""))
# [1] "Number of checked genes: 13"
# > print(paste("Number of checked genes after removing genes with low expression: ", nrow(counts[which(counts$Geneid %in% geneList),]), sep = ""))
# [1] "Number of checked genes after removing genes with low expression: 10"
# > if(nrow(counts[which(counts$Geneid %in% geneList),]) < length(geneList)) {
#   +     print("Genes are removed due to low expression:")
#   +     print(geneList[which(!(geneList %in% counts$Geneid))])
#   + } else {
#     +     print("All genes are kept after removing genes with low expression")
#     + }
# [1] "Genes are removed due to low expression:"
# [1] "ALB"     "S100A12" "SNCA"   

geneList <- geneList[which(geneList %in% counts$Geneid)]

group<-colnames(counts)[2:ncol(counts)]
group
# > group
# [1] "MB231_CON_1_S17_L001.sorted.bam"    "MB231_CON_2_S7_L001.sorted.bam"     "MB231_CON_3_S13_L001.sorted.bam"    "MB231_TEPA24_1_S8_L001.sorted.bam" 
# [5] "MB231_TEPA24_2_S4_L001.sorted.bam"  "MB231_TEPA24_3_S11_L001.sorted.bam" "MB231_TEPA8_1_S2_L001.sorted.bam"   "MB231_TEPA8_2_S10_L001.sorted.bam" 
# [9] "MB231_TEPA8_3_S14_L001.sorted.bam" 
# group<-gsub(pattern="[.][0-9]$",replacement="",group)
group<-gsub(pattern="[_][0-9][_].*",replacement="",group)
group
# > group
# [1] "MB231_CON"    "MB231_CON"    "MB231_CON"    "MB231_TEPA24" "MB231_TEPA24" "MB231_TEPA24" "MB231_TEPA8"  "MB231_TEPA8"  "MB231_TEPA8"

# Create colors for groups
colGrp<-group
for(i in 1:length(colGrp)){
  if(colGrp[i] == "MB231_CON"){
    colGrp[i] <- "goldenrod"
  }
  if(colGrp[i] == "MB231_TEPA24"){
    colGrp[i] <- "darkgreen"
  }
  if(colGrp[i] == "MB231_TEPA8"){
    colGrp[i] <- "lightskyblue1"
  }
}

#================================================================
# PCA analysis
#================================================================
# tempTpm<-tpm[,2:ncol(tpm)]
# tempTpm[tempTpm==0]<-0.00001
# tTpm<-t(tempTpm)
# logTpm<-log(tTpm)
# 
# tpm.pca<-prcomp(logTpm)
# dataPCA<-cbind(logTpm,grouping=group)
# autoplot(tpm.pca,data=dataPCA,colour='grouping',label=T)
# filename<-paste(outDir, "/RNAseq_PCA",".png",sep="")
# png(filename,res=100,height=800,width=1000)
# autoplot(tpm.pca,data=dataPCA,colour='grouping')
# dev.off()

y<-DGEList(counts=counts[,2:ncol(counts)],group=group,genes=counts[,1])
y<-calcNormFactors(y)

pdf(file=paste(outDir, "/QC_metrics_plots.pdf", sep=""))
par(mar=c(5.1,4.1,4.1,10))
plotMDS(y,main="MDS plot of RNA-seq data")
plotMD(y,column=1)
y<-estimateCommonDisp(y)
y<-estimateGLMTrendedDisp(y)
y<-estimateTagwiseDisp(y)
plotBCV(y)
dev.off()

#================================================================
# Perform DEGs analysis
#================================================================
design<-model.matrix(~ -1+factor(group))
colnames(design) <- levels(factor(group))

# TEPA24
contrast.matrix<-makeContrasts(DrugvsControl.MB231.TEPA24=MB231_TEPA24-MB231_CON,levels=design)
# TEPA8
contrast.matrix<-cbind(contrast.matrix, makeContrasts(DrugvsControl.MB231.TEPA8=MB231_TEPA8-MB231_CON,levels=design))

fit<-glmFit(y,design)

for(i in 1:ncol(contrast.matrix)){
  deg_lrt<-glmLRT(fit,contrast=contrast.matrix[,i])
  is.de<-decideTestsDGE(deg_lrt)
  plotMD(deg_lrt,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
  deg_all<-topTags(deg_lrt,n=nrow(deg_lrt),adjust.method="BH",sort.by="PValue")
  deg_all_table<-as.data.frame(deg_all)
  FC<-logratio2foldchange(deg_all_table$logFC,base=2)
  deg_expTab<-cbind.data.frame(deg_all_table,FC)
  m<-match(deg_expTab$genes,egSymbol$symbol)
  deg_expTab$entrez_id<-egSymbol$gene_id[m]
  
  # Only get interested genes
  deg_expTab <- deg_expTab[which(deg_expTab$genes %in% geneList),]
  
  #Add mygene summaries
  geneSummaries<-as.data.frame(getGenes(geneids=c(deg_expTab$entrez_id),fields="all",return.as="DataFrame")[,c("query","summary")])
  m2<-match(deg_expTab$entrez_id,geneSummaries$query)
  deg_expTab$summary<-geneSummaries$summary[m2]
  
  fileName<-paste(outDir, "/", colnames(contrast.matrix)[i],".txt",sep="")
  write.table(deg_expTab,fileName,sep="\t",row.names=F)
  # temp<-merge(deg_expTab,tpm,by="genes")
  # reOrdCol<-c(1,8,2:7,9:ncol(temp))
  # deg_expTab_tpm<-temp[order(temp$PValue),reOrdCol]
  # fileName<-paste(outDir, "/", colnames(contrast.matrix)[i],"_tpm.txt",sep="")
  # write.table(deg_expTab_tpm,fileName,sep="\t",row.names=F)
  
  #Enrichment analysis with clusterProfiler/DAVID API
  upRegGenes<-deg_expTab$entrez_id[which(deg_expTab$FC>=FC_threshold&deg_expTab$PValue<=0.05)]
  downRegGenes<-deg_expTab$entrez_id[which(deg_expTab$FC<=-FC_threshold&deg_expTab$PValue<=0.05)]
  print(deg_expTab$genes[which(deg_expTab$FC>=FC_threshold&deg_expTab$PValue<=0.05)])
  print(deg_expTab$genes[which(deg_expTab$FC<=-FC_threshold&deg_expTab$PValue<=0.05)])
  
  #Volcano Plot with points FC>=|2| AND FDR<0.05 highlighted
  fileName<-paste(outDir, "/", colnames(contrast.matrix)[i],"_volcanoPlot.pdf",sep="")
  pdf(fileName,  width = 6, height = 6, onefile=FALSE)
  title<-paste(colnames(contrast.matrix)[i],": FC vs FDR",sep="")
  with(deg_expTab,plot(deg_expTab$logFC,-log10(deg_expTab$FDR),pch=20,main=title, xlab="log2FoldChange",ylab="-log10(FDR)"))
  legend("topright", legend=c("FDR < 0.05", "FDR < 0.05 and absolute(logFC) > 1)", "FDR >= 0.05 or absolute(logFC) <= 1)"), col = c("red", "green", "black"), pch=20)
  with(deg_expTab, text(deg_expTab$logFC,-log10(deg_expTab$FDR), labels=genes, cex= 0.9))
  with(subset(deg_expTab,FDR<0.05),points(logFC,-log10(FDR),pch=20,col="red"))
  with(subset(deg_expTab,abs(logFC)>1),points(logFC,-log10(FDR),pch=20,col="blue"))
  with(subset(deg_expTab,FDR<0.05&abs(logFC)>1),points(logFC,-log10(FDR),pch=20,col="green"))
  dev.off()
  
  # # Filter biologically significant results FC>=|2|
  # deg_FC2 <- deg_expTab[which(abs(deg_expTab$FC) >= 2),]
  # deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
  # if(nrow(deg_FC2) > 1){
  #   
  #   gene_list <- deg_FC2$genes
  #   mat <- as.matrix(tpm[,2:ncol(tpm)])
  #   rownames(mat)<-tpm$genes
  #   mat[mat==0] <- 0.001
  #   gn<-rownames(mat)
  #   index<-which(gn %in% gene_list)
  #   tpm_FC2<-as.matrix(mat[index,])
  #   
  #   fileName<-paste(outDir, "/", colnames(contrast.matrix)[i],"_SigGenes_FC2.png",sep="")
  #   title<-paste(colnames(contrast.matrix)[i]," FC>=|2|",sep="")
  #   png(fileName, height=700, width=800,res=100)
  #   heatmap.2(x=log2(tpm_FC2),col=bluered(256),scale="row",trace="none",
  #             cexCol=0.9,cexRow=0.7,srtCol=35,main=title,
  #             labRow=rownames(tpm_FC2),colCol=colGrp,ColSideColors=colGrp)
  #   dev.off()
  # }
  # # Filter biologically significant results FC>=|2| and statistically significant results fdr<0.05
  # deg_FC2 <- deg_expTab[which((abs(deg_expTab$FC)>=2) & (deg_expTab$FDR<0.05)),]
  # deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
  # if(nrow(deg_FC2) > 1){
  #   gene_list <- deg_FC2$genes
  #   mat <- as.matrix(tpm[,2:ncol(tpm)])
  #   rownames(mat)<-tpm$genes
  #   mat[mat==0] <- 0.001
  #   gn<-rownames(mat)
  #   index<-which(gn %in% gene_list)
  #   tpm_FC2<-as.matrix(mat[index,])
  #   
  #   fileName<-paste(outDir, "/", colnames(contrast.matrix)[i],"_SigGenes_FC2_sigFDR.png",sep="")
  #   title<-paste(colnames(contrast.matrix)[i]," FC>=|2| & FDR<0.05",sep="")
  #   png(fileName,width=800,height=700,res=100)
  #   heatmap.2(x=log2(tpm_FC2),col=bluered(256),scale="row",trace="none",
  #             cexCol=0.9,cexRow=0.7,srtCol=35,main=title,
  #             labRow=rownames(tpm_FC2),colCol=colGrp,ColSideColors=colGrp)
  #   dev.off()
  # }
  
  # # Heatmap for all genes
  # deg_FC2 <- deg_expTab
  # deg_FC2 <- deg_FC2[order(deg_FC2$FC, decreasing = FALSE),]
  # if(nrow(deg_FC2) > 1){
  #   gene_list <- deg_FC2$genes
  #   mat <- as.matrix(tpm[,2:ncol(tpm)])
  #   rownames(mat)<-tpm$genes
  #   mat[mat==0] <- 0.001
  #   gn<-rownames(mat)
  #   index<-which(gn %in% gene_list)
  #   tpm_FC2<-as.matrix(mat[index,])
  #   
  #   fileName<-paste(outDir, "/", "AllCheckedGenes.png",sep="")
  #   title<-paste("All checked genes",sep="")
  #   png(fileName,width=800,height=700,res=100)
  #   heatmap.2(x=log2(tpm_FC2),col=bluered(256),scale="row",trace="none",
  #             cexCol=0.9,cexRow=0.7,srtCol=35,main=title,
  #             labRow=rownames(tpm_FC2),colCol=colGrp,ColSideColors=colGrp)
  #   dev.off()
  # }
  # 
}
