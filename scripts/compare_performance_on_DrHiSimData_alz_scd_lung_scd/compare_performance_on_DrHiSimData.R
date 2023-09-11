#! /usr/bin/Rscript

# compare_performance_on_DrHiSimData.R
# Juok Cho
# 2021 04 19



compare_performance_on_DrHiSimData <- function(DataN, filter_rate){
  
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(sctransform)
  library(Linnorm)
  library(ggplot2)
  library("SingleCellExperiment")
  
  source("./scripts/compare_performance_on_DrHiSimData/comp_ARI.ASW.N2S.UMAP.DEGR_DrHiSimData.R")
  
  file.path.result = "./results/compare_performance_on_DrHiSimData/"
  
  ############# load real data ######################
  n = as.numeric(DataN)
  filter_rate = filter_rate
  
  #DrHiDepthSimData v3&v4
  path.data = "./data/DrHiSimData/"
  file.list = dir(path.data) 
  w = grep(".data.", file.list)
  datafile.list = file.list[w]
  dataName.list = gsub(".data.RDS", "", datafile.list)
  
  # raw count Data 
  #class: SingleCellExperiment
  fileName = paste(path.data, datafile.list[n], sep="")
  simData <- readRDS(fileName) 
  count.data = counts(simData)
  real.subtype <- simData$Group   
  rownames(count.data)=gsub("_","-", rownames(count.data))

  # degData
  degfileName = paste(path.data, dataName.list[n], ".degs.RDS", sep="")
  degData <- readRDS(degfileName) 
  
  message(paste("n:",n ," - ", dataName.list[n], " ...load realData"))
  message(paste("filter_rate:", filter_rate))
  
  # all zero genes filtering
  rowsum = rowSums(count.data)
  w.rm = which(rowsum==0)
  if(length(w.rm)>0){ count.data = count.data[-w.rm, ]}


  # generate ARI table
  results = comp_ARI.realData(dataName=dataName.list[n], count.data = count.data, real.subtype = real.subtype, filter_rate, degData=degData)
  
  ARI.table = as.matrix(results$ARI.mat)
  ASW.table = as.matrix(results$ASW.mat)
  DEGall.table = as.matrix(results$DEG.mat.all)
  DEGup.table = as.matrix(results$DEG.mat.up)
  DEGdown.table = as.matrix(results$DEG.mat.down)
  EGR.G1.table = as.matrix(results$EGR.mat.G1)
  EGR.G2.table = as.matrix(results$EGR.mat.G2)
  EGR.G3.table = as.matrix(results$EGR.mat.G3)
  EGR.G4.table = as.matrix(results$EGR.mat.G4)
  EGR.G5.table = as.matrix(results$EGR.mat.G5)

  tableName = paste("DrHiDepth",dataName.list[n], filter_rate, sep="_")

  write.table(ARI.table,
              paste(file.path.result, "ARI.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
  
  write.table(ASW.table,
              paste(file.path.result, "ASW.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(DEGall.table,
              paste(file.path.result, "DEGall.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(DEGup.table,
              paste(file.path.result, "DEGup.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(DEGdown.table,
              paste(file.path.result, "DEGdown.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(EGR.G1.table,
              paste(file.path.result, "EGR.G1.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(EGR.G2.table,
              paste(file.path.result, "EGR.G2.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(EGR.G3.table,
              paste(file.path.result, "EGR.G3.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(EGR.G4.table,
              paste(file.path.result, "EGR.G4.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

  write.table(EGR.G5.table,
              paste(file.path.result, "EGR.G5.table", tableName,".txt", sep=""),
              row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
 
  message(paste("saved the ARI,ASW,DEG,EGR table ..."))
  
}


###### invoke r function when job is submitted ##########
Args <- commandArgs(TRUE);
compare_realData(DataN=Args[1], filter_rate=Args[2])