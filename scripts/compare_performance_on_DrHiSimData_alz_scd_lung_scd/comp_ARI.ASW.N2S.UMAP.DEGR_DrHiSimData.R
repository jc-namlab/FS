# comp_ARI.ASW.N2S.UMAP.DEGR_DrHiSimData.R
# Juok Cho
# 2021 11 01


comp_ARI.ASW.N2S.UMAP.DEGR_DrHiSimData = function(dataName, count.data, real.subtype, filter_rate, degData=degData){
  nfeatures.list = c(500, 1000, 2000, 4000)
  tmp = matrix(0, ncol=22, nrow=4)
  colnames(tmp) =   c("scTv2.HVG.vst","scTv2.HDG", "scTv2.HEG.cpm",
                      "scTv2.HVG.plus.HDG", "scTv2.HVG.plus.HEG",
                      "scTv2.NBDrop", "scTv2.NBDisp", "scTv2.M3Drop", 
                      "scTv2.HIPPO","scTv2.HVG.disp","scTv2.DUBStepR",
                      "logN.HVG.vst", "logN.HDG", "logN.HEG.cpm",
                      "logN.HVG.plus.HDG", "logN.HVG.plus.HEG",
                      "logN.NBDrop", "logN.NBDisp", "logN.M3Drop",
                      "logN.HIPPO",  "logN.HVG.disp", "logN.DUBStepR")
  
                     
  
  rownames(tmp) = nfeatures.list
  
  ARI.mat = tmp
  ASW.mat = tmp
  DEG.mat.all = tmp
  DEG.mat.up = tmp
  DEG.mat.down = tmp
  EGR.mat.G1 = tmp
  EGR.mat.G2 = tmp
  EGR.mat.G3 = tmp
  EGR.mat.G4 = tmp
  EGR.mat.G5 = tmp
  
  # cpm exp level 
  cpm.data = edgeR::cpm(count.data)
  gmean = rowMeans(cpm.data)
  q = quantile(gmean, c(0.2,0.4,0.6,0.8))
  w1 = which(gmean < q[1])
  w2 = which(gmean <q[2] & gmean >=q[1])
  w3 = which(gmean < q[3] & gmean >=q[2])
  w4 = which(gmean <q[4] & gmean >=q[3] )
  w5 = which(gmean >= q[4])
  expGroups = list()
  expGroups$G1 = names(gmean)[w1]
  expGroups$G2 = names(gmean)[w2]
  expGroups$G3 = names(gmean)[w3]
  expGroups$G4 = names(gmean)[w4]
  expGroups$G5 = names(gmean)[w5]
  
  # filter_rate: filter low 10 % , filter low 20%
  # cpm.data = edgeR::cpm(count.data)
  # gmean = rowMeans(cpm.data)
  cutoff = quantile(gmean, as.numeric(filter_rate))
  w <- which(gmean >= cutoff)
  count.data <- count.data[w, ]
  

  # create seurat obj
  N.minCells = 0
  obj.org <- CreateSeuratObject(counts = count.data, project = "real", min.cells = N.minCells)
  obj.org[["percent.mt"]] <- PercentageFeatureSet(obj.org, pattern = "^MT-")
  message("created seurat obj...")
  obj.org
  
  obj = obj.org
  
  # # QC & filter cells
  # obj <- subset(obj.org, subset = nFeature_RNA > 200 & percent.mt < 5)
  # message("obj after QC & filter cells - ^MT-...")
  # obj
  
  # n: nfeatures.list = c(500, 1000, 2000, 4000)
  N.f = length(nfeatures.list)
  for(k in 1:N.f){
    print(paste("k:",k,sep=""))
    
    k.results = comp_N.features(dataName, nfeatures=nfeatures.list[k], obj, real.subtype, expGroups=expGroups, degData=degData)
    
    ARI.mat[k,] =k.results$ARI.mat
    ASW.mat[k,] =k.results$ASW.mat
    
    DEG.mat.all[k,] = k.results$DEG.mat.all
    DEG.mat.up[k,] = k.results$DEG.mat.up
    DEG.mat.down[k,] = k.results$DEG.mat.down
    
    EGR.mat.G1[k,] = k.results$EGR.mat.G1
    EGR.mat.G2[k,] = k.results$EGR.mat.G2
    EGR.mat.G3[k,] = k.results$EGR.mat.G3
    EGR.mat.G4[k,] = k.results$EGR.mat.G4
    EGR.mat.G5[k,] = k.results$EGR.mat.G5

  }
  
  results = list()
  results$ARI.mat = ARI.mat
  results$ASW.mat = ASW.mat
  results$DEG.mat.all = DEG.mat.all
  results$DEG.mat.up = DEG.mat.up
  results$DEG.mat.down = DEG.mat.down
  results$EGR.mat.G1 = EGR.mat.G1
  results$EGR.mat.G2 = EGR.mat.G2
  results$EGR.mat.G3 = EGR.mat.G3
  results$EGR.mat.G4 = EGR.mat.G4
  results$EGR.mat.G5 = EGR.mat.G5
  
  return(results)

}




comp_N.features = function(dataName, nfeatures, obj, real.subtype, expGroups=expGroups, degData=degData){
  
  library(aricode)
  library(Seurat)
  library(sctransform)
  library(Linnorm)
  library(edgeR)
  library(DESeq2)
  library(scran)
  library(cluster)
  library("glmGamPoi")
  library("sctransformV2Test")
  library(M3Drop)
  library(HIPPO)
  library(DUBStepR.TEST)
  
  print(paste("nfeatures:", nfeatures, sep=""))
  
  tmp = matrix(0, ncol=22, nrow=1)
  colnames(tmp) =   c("scTv2.HVG.vst","scTv2.HDG", "scTv2.HEG.cpm",
                      "scTv2.HVG.plus.HDG", "scTv2.HVG.plus.HEG",
                      "scTv2.NBDrop", "scTv2.NBDisp", "scTv2.M3Drop", 
                      "scTv2.HIPPO","scTv2.HVG.disp","scTv2.DUBStepR",
                      "logN.HVG.vst", "logN.HDG", "logN.HEG.cpm",
                      "logN.HVG.plus.HDG", "logN.HVG.plus.HEG",
                      "logN.NBDrop", "logN.NBDisp", "logN.M3Drop",
                      "logN.HIPPO",  "logN.HVG.disp", "logN.DUBStepR")
  
  
  ARI.mat = tmp
  ASW.mat = tmp
  
  DEG.mat.all = tmp
  DEG.mat.up = tmp
  DEG.mat.down = tmp
  
  EGR.mat.G1 = tmp
  EGR.mat.G2 = tmp
  EGR.mat.G3 = tmp
  EGR.mat.G4 = tmp
  EGR.mat.G5 = tmp
  
  #***** Run different normalization method
  
  #1.sctransform::vst - scT.v2 HVG.counts
  obj.scTv2.HVG.vst = obj
  scT.v2_out <- sctransformV2Test::vst_S2N(obj[["RNA"]]@counts, vst.flavor = "v2")
  obj.scTv2.HVG.vst[["RNA"]]@data <- scT.v2_out$y
  
  genesToUse = rownames(scT.v2_out$y)
  cellsToUse = colnames(scT.v2_out$y)
  obj[["RNA"]]@counts = obj[["RNA"]]@counts[genesToUse, cellsToUse]
  obj[["RNA"]]@data = obj[["RNA"]]@data[genesToUse, cellsToUse]
  obj.scTv2.HVG.vst[["RNA"]]@counts= obj[["RNA"]]@counts
  
  obj.scTv2.HVG.vst <- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.HVG.vst, subtype = real.subtype, method="HVG.vst",  nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.HVG.vst))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.HVG.vst))
  i=1; j=1
  ARI.mat[i,j] = obj.scTv2.HVG.vst@misc$ARI
  ASW.mat[i,j] = obj.scTv2.HVG.vst@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method1.scT v2.HVG.vst counts.. ARI Done!")

  
  #2. scT v2.HDG
  obj.scTv2.HDG = obj.scTv2.HVG.vst
  obj.scTv2.HDG <- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.HDG, subtype = real.subtype,  method="HDG", nfeatures = nfeatures, dist="poisson")
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.HDG))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.HDG))
  i=1; j=2
  ARI.mat[i,j] = obj.scTv2.HDG@misc$ARI
  ASW.mat[i,j] = obj.scTv2.HDG@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method2.scT v1.HDG...ARI Done!")

  #3. scTv2.HEG.cpm
  obj.scTv2.HEG.cpm = obj.scTv2.HVG.vst
  obj.scTv2.HEG.cpm<- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.HEG.cpm, subtype = real.subtype, method="HEG.cpm", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.HEG.cpm))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.HEG.cpm))
  i=1; j=3
  ARI.mat[i,j] = obj.scTv2.HEG.cpm@misc$ARI
  ASW.mat[i,j] = obj.scTv2.HEG.cpm@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method3.scT v1.HEG.cpm...ARI Done!")
  
  
  #4. "scTv2.HVG.plus.HDG"
  obj.scTv2.HVG.plus.HDG = obj.scTv2.HVG.vst
  obj.scTv2.HVG.plus.HDG<- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.HVG.plus.HDG, subtype = real.subtype, method="HVG.plus.HDG", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.HVG.plus.HDG))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.HVG.plus.HDG))
  i=1; j=4
  ARI.mat[i,j] = obj.scTv2.HVG.plus.HDG@misc$ARI
  ASW.mat[i,j] = obj.scTv2.HVG.plus.HDG@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method4.scTv2.HVG.plus.HDG...ARI Done!")
  
  #5. "scTv2.HVG.plus.HEG"
  obj.scTv2.HVG.plus.HEG = obj.scTv2.HVG.vst
  obj.scTv2.HVG.plus.HEG<- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.HVG.plus.HEG, subtype = real.subtype, method="HVG.plus.HEG", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.HVG.plus.HEG))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.HVG.plus.HEG))
  i=1; j=5
  ARI.mat[i,j] = obj.scTv2.HVG.plus.HEG@misc$ARI
  ASW.mat[i,j] = obj.scTv2.HVG.plus.HEG@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method5.scTv2.HVG.plus.HEG...ARI Done!")
  
  #6. "scTv2.NBDrop"
  obj.scTv2.NBDrop = obj.scTv2.HVG.vst
  obj.scTv2.NBDrop<- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.NBDrop, subtype = real.subtype, method="NBDrop", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.NBDrop))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.NBDrop))
  i=1; j=6
  ARI.mat[i,j] = obj.scTv2.NBDrop@misc$ARI
  ASW.mat[i,j] = obj.scTv2.NBDrop@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method6.scTv2.NBDrop...ARI Done!")
  
  #7. "scTv2.NBDisp"
  obj.scTv2.NBDisp = obj.scTv2.HVG.vst
  obj.scTv2.NBDisp<- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.NBDisp, subtype = real.subtype, method="NBDisp", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.NBDisp))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.NBDisp))
  i=1; j=7
  ARI.mat[i,j] = obj.scTv2.NBDisp@misc$ARI
  ASW.mat[i,j] = obj.scTv2.NBDisp@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method7.scTv2.NBDisp...ARI Done!")
  
  #8. "scTv2.M3Drop"
  obj.scTv2.M3Drop = obj.scTv2.HVG.vst
  obj.scTv2.M3Drop<- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.M3Drop, subtype = real.subtype, method="M3Drop", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.M3Drop))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.M3Drop))
  i=1; j=8
  ARI.mat[i,j] = obj.scTv2.M3Drop@misc$ARI
  ASW.mat[i,j] = obj.scTv2.M3Drop@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method8.scTv2.M3Drop...ARI Done!")
  
  #9. "scTv2.HIPPO"
  obj.scTv2.HIPPO = obj.scTv2.HVG.vst
  obj.scTv2.HIPPO<- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.HIPPO, subtype = real.subtype, method="HIPPO", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.HIPPO))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.HIPPO))
  i=1; j=9
  ARI.mat[i,j] = obj.scTv2.HIPPO@misc$ARI
  ASW.mat[i,j] = obj.scTv2.HIPPO@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method9.scTv2.HIPPO...ARI Done!")
  
  #10. scT v2.HVG.disp
  obj.scTv2.HVG.disp = obj.scTv2.HVG.vst
  obj.scTv2.HVG.disp <- generate_clustering_results.Vresol.obj_features(obj.name = obj.scTv2.HVG.disp, subtype = real.subtype,  method="HVG.disp", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.HVG.disp))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.HVG.disp))
  i=1; j=10
  ARI.mat[i,j] = obj.scTv2.HVG.disp@misc$ARI
  ASW.mat[i,j] = obj.scTv2.HVG.disp@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method10.scT v2 HVG.disp...ARI Done!")
  
  #11. scT v2 DUBStepR
  obj.scTv2.DUBStepR = obj.scTv2.HVG.vst
  #obj.scTv2.DUBStepR <- generate_clustering_results.Vresol.obj_features(model.name = "scTv2",obj.name = obj.scTv2.DUBStepR, subtype = real.subtype,  method="DUBStepR", nfeatures = nfeatures, dist=NULL)
  obj.scTv2.DUBStepR <- generate_clustering_results.Vresol.obj_features_DUBStepR(model.name = "scTv2",dataName, obj.name = obj.scTv2.DUBStepR, subtype = real.subtype,  method="DUBStepR", nfeatures = nfeatures, dist=NULL, filterN="filter0.05")
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.scTv2.DUBStepR))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.scTv2.DUBStepR))
  i=1; j=11
  ARI.mat[i,j] = obj.scTv2.DUBStepR@misc$ARI
  ASW.mat[i,j] = obj.scTv2.DUBStepR@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method11.scT v2 DUBStepR...ARI Done!")
  

  
  
  # 12. logN.HVG.vst
  obj.logN.HVG.vst = obj
  obj.logN.HVG.vst <- NormalizeData(obj.logN.HVG.vst, normalization.method = "LogNormalize", scale.factor = 10000)
  obj.logN.HVG.vst <- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.HVG.vst, subtype = real.subtype, method="HVG.vst", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.HVG.vst))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.HVG.vst))
  i=1; j=12
  ARI.mat[i,j] = obj.logN.HVG.vst@misc$ARI
  ASW.mat[i,j] = obj.logN.HVG.vst@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method12.logN.HVG.vst...ARI Done!")
  

  #13. logN.HDG
  obj.logN.HDG = obj.logN.HVG.vst
  obj.logN.HDG<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.HDG, subtype = real.subtype, method="HDG", nfeatures = nfeatures, dist="poisson")
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.HDG))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.HDG))
  i=1; j=13
  ARI.mat[i,j] = obj.logN.HDG@misc$ARI
  ASW.mat[i,j] = obj.logN.HDG@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method13.logN.HDG...ARI Done!")
  
  
  #14.logN.HEG.cpm
  obj.logN.HEG.cpm = obj.logN.HVG.vst
  obj.logN.HEG.cpm<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.HEG.cpm, subtype = real.subtype, method="HEG.cpm", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.HEG.cpm))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.HEG.cpm))
  i=1; j=14
  ARI.mat[i,j] = obj.logN.HEG.cpm@misc$ARI
  ASW.mat[i,j] = obj.logN.HEG.cpm@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method14.logN.HEG.cpm...ARI Done!")
  
  #15. "logN.HVG.plus.HDG"
  obj.logN.HVG.plus.HDG = obj.logN.HVG.vst
  obj.logN.HVG.plus.HDG<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.HVG.plus.HDG, subtype = real.subtype, method="HVG.plus.HDG", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.HVG.plus.HDG))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.HVG.plus.HDG))
  i=1; j=15
  ARI.mat[i,j] = obj.logN.HVG.plus.HDG@misc$ARI
  ASW.mat[i,j] = obj.logN.HVG.plus.HDG@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method15.logN.HVG.plus.HDG...ARI Done!")
  
  #16. "logN.HVG.plus.HEG"
  obj.logN.HVG.plus.HEG = obj.logN.HVG.vst
  obj.logN.HVG.plus.HEG<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.HVG.plus.HEG, subtype = real.subtype, method="HVG.plus.HEG", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.HVG.plus.HEG))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.HVG.plus.HEG))
  i=1; j=16
  ARI.mat[i,j] = obj.logN.HVG.plus.HEG@misc$ARI
  ASW.mat[i,j] = obj.logN.HVG.plus.HEG@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method16.logN.HVG.plus.HEG...ARI Done!")
  
  #17. "logN.NBDrop"
  obj.logN.NBDrop = obj.logN.HVG.vst
  obj.logN.NBDrop<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.NBDrop, subtype = real.subtype, method="NBDrop", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.NBDrop))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.NBDrop))
  i=1; j=17
  ARI.mat[i,j] = obj.logN.NBDrop@misc$ARI
  ASW.mat[i,j] = obj.logN.NBDrop@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method17.logN.NBDrop...ARI Done!")
  
  #18. "logN.NBDisp"
  obj.logN.NBDisp = obj.logN.HVG.vst
  obj.logN.NBDisp<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.NBDisp, subtype = real.subtype, method="NBDisp", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.NBDisp))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.NBDisp))
  i=1; j=18
  ARI.mat[i,j] = obj.logN.NBDisp@misc$ARI
  ASW.mat[i,j] = obj.logN.NBDisp@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method18.logN.NBDisp...ARI Done!")
  
  #19. "logN.M3Drop"
  obj.logN.M3Drop = obj.logN.HVG.vst
  obj.logN.M3Drop<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.M3Drop, subtype = real.subtype, method="M3Drop", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.M3Drop))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.M3Drop))
  i=1; j=19
  ARI.mat[i,j] = obj.logN.M3Drop@misc$ARI
  ASW.mat[i,j] = obj.logN.M3Drop@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method19.logN.M3Drop...ARI Done!")
  
  
  #20. "logN.HIPPO"
  obj.logN.HIPPO = obj.logN.HVG.vst
  obj.logN.HIPPO<- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.HIPPO, subtype = real.subtype, method="HIPPO", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.HIPPO))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.HIPPO))
  i=1; j=20
  ARI.mat[i,j] = obj.logN.HIPPO@misc$ARI
  ASW.mat[i,j] = obj.logN.HIPPO@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method20.logN.HIPPO...ARI Done!")
  
  #21. logN.HVG.disp
  obj.logN.HVG.disp = obj.logN.HVG.vst
  obj.logN.HVG.disp <- generate_clustering_results.Vresol.obj_features(obj.name = obj.logN.HVG.disp, subtype = real.subtype,  method="HVG.disp", nfeatures = nfeatures, dist=NULL)
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.HVG.disp))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.HVG.disp))
  i=1; j=21
  ARI.mat[i,j] = obj.logN.HVG.disp@misc$ARI
  ASW.mat[i,j] = obj.logN.HVG.disp@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method21.logN HVG.disp...ARI Done!")
  
 
  #22. logN.DUBStepR
  obj.logN.DUBStepR = obj.logN.HVG.vst
  obj.logN.DUBStepR <- generate_clustering_results.Vresol.obj_features_DUBStepR(model.name = "logN",dataName,  obj.name = obj.logN.DUBStepR, subtype = real.subtype,  method="DUBStepR", nfeatures = nfeatures, dist=NULL,filterN="filter0.05")
  EGR <- featuresExpGroupsV2(expGroups, features = VariableFeatures(obj.logN.DUBStepR))
  DEGR <- trueDEGRV2(degData, features= VariableFeatures(obj.logN.DUBStepR))
  i=1; j=22
  ARI.mat[i,j] = obj.logN.DUBStepR@misc$ARI
  ASW.mat[i,j] = obj.logN.DUBStepR@misc$ASW
  DEG.mat.all[i,j] = DEGR$all
  DEG.mat.up[i,j] = DEGR$up
  DEG.mat.down[i,j] = DEGR$down
  EGR.mat.G1[i,j] = EGR$G1
  EGR.mat.G2[i,j] = EGR$G2
  EGR.mat.G3[i,j] = EGR$G3
  EGR.mat.G4[i,j] = EGR$G4
  EGR.mat.G5[i,j] = EGR$G5
  print("method22.logN DUBStepR...ARI Done!")
 
 
  
  #######################################
  
  message("finish! return ARI.mat")
  
  results = list()
  results$ARI.mat = ARI.mat
  results$ASW.mat = ASW.mat
  results$DEG.mat.all = DEG.mat.all
  results$DEG.mat.up = DEG.mat.up
  results$DEG.mat.down = DEG.mat.down
  results$EGR.mat.G1 = EGR.mat.G1
  results$EGR.mat.G2 = EGR.mat.G2
  results$EGR.mat.G3 = EGR.mat.G3
  results$EGR.mat.G4 = EGR.mat.G4
  results$EGR.mat.G5 = EGR.mat.G5

  return(results)
}



generate_clustering_results.basic = function(obj.name, ari.clus.resol=0.7){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  
  #***** 4. Identification of highly variable features (feature selection)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
  #***** 5. Scaling the data 
  # (,so that highly-expressed genes do not dominate)
  
  tryCatch({
    test <- ScaleData(obj, features = rownames(obj), verbose =FALSE)},
    error = function(e) print("ScaleData - initial trial returns error...Let's try again!"),
    warning = function(w) print("ScaleData - initial trial returns warning...Let's try again!"),
    finally = NULL)
  
  obj <- ScaleData(obj, features = rownames(obj), verbose =FALSE)
  
  #***** 6. Perform linear dimensional reduction (PCA)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE)
  
  #***** 7. Cluster the cells
  obj <- FindNeighbors(obj, dims = 1:10, k.param=20, n.trees=50, verbose =FALSE) ## (Jaccard similarity)
  obj <- FindClusters(obj, resolution = ari.clus.resol, verbose =FALSE) ##previously, resolution = 0.5, apply modularity optimization techniques to iteratively group cells together 
  
  return(obj$seurat_clusters)
}

generate_clustering_results.Vresol <- function(obj.name, subtype){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(subtype))
  #***** 4. Identification of highly variable features (feature selection)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose =FALSE)
  #***** 5. Scaling the data 
  # (,so that highly-expressed genes do not dominate)
  
  tryCatch({
    test <- ScaleData(obj, features = rownames(obj), verbose =FALSE)},
    error = function(e) print("ScaleData - initial trial returns error...Let's try again!"),
    warning = function(w) print("ScaleData - initial trial returns warning...Let's try again!"),
    finally = NULL)
  
  obj <- ScaleData(obj, features = rownames(obj), verbose =FALSE)
  
  #***** 6. Perform linear dimensional reduction (PCA)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE)
  
  #***** 7. Cluster the cells
  obj <- FindNeighbors(obj, dims = 1:10, k.param=20, n.trees=50, verbose =FALSE) ## (Jaccard similarity)
  
  ari.clus.resol.list = c(seq(1.2,0.1, -0.2), 0.1)
  N = length(ari.clus.resol.list)
  message(paste("ARI.clus.update:TRUE --- try different resolution values... in [1.2 - 0.1 ]"))
  
  for(i in 1:N){
    tmp <- FindClusters(obj, resolution = ari.clus.resol.list[i], verbose =FALSE) ##previously, resolution = 0.5, apply modularity optimization techniques to iteratively group cells together 
    clus.update = tmp$seurat_clusters
    N.clus <- length(table(clus.update))
    message(paste("i:",i, "N.clus:", N.clus, "clus.minNum:", clus.minNum))
    
    if(N.clus > clus.minNum){ 
      N.clus.prev <- N.clus
      clus.update.prev <- clus.update
    }else if(N.clus == clus.minNum){
      message(paste("i:", i, "ari.clus.resol:", ari.clus.resol.list[i], "N.clus:", N.clus, "clus.minNum:", clus.minNum ))
      break;
      
    }else if(N.clus < clus.minNum){
      if(i!=1){
        N.clus.next <- N.clus
        list.N.clus <- c(N.clus.prev, N.clus.next)
        w <- which.min(abs(list.N.clus - clus.minNum))
        if(w==1){clus.update <- clus.update.prev}
      }
      break;
    }
  }
  return(clus.update)
}



generate_clustering_results.obj <- function(obj.name, subtype, fixed.resol=fixed.resol.value){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(subtype))
  #***** 4. Identification of highly variable features (feature selection)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose =FALSE)
  #***** 5. Scaling the data 
  # (,so that highly-expressed genes do not dominate)
  
  tryCatch({
    test <- ScaleData(obj, features = rownames(obj), verbose =FALSE)},
    error = function(e) print("ScaleData - initial trial returns error...Let's try again!"),
    warning = function(w) print("ScaleData - initial trial returns warning...Let's try again!"),
    finally = NULL)
  
  obj <- ScaleData(obj, features = rownames(obj), verbose =FALSE)
  
  #***** 6. Perform linear dimensional reduction (PCA)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE)
  
  #***** 7. Cluster the cells
  N.PCs <- length(obj[['pca']])
  obj <- FindNeighbors(obj, dims = 1:N.PCs, k.param=20, n.trees=50, verbose =FALSE) ## (Jaccard similarity)
  
  ari.clus.resol.list = fixed.resol #c(1.2)
  N = length(ari.clus.resol.list)
  message(paste("ARI.clus.update:TRUE --- try different resolution values... in [1.2 -0.1, 0.01 ]"))
  
  for(i in 1:N){
    tmp <- FindClusters(obj, resolution = ari.clus.resol.list[i], verbose =FALSE) ##previously, resolution = 0.5, apply modularity optimization techniques to iteratively group cells together 
    clus.update = tmp$seurat_clusters
    N.clus <- length(table(clus.update))
    message(paste("i:",i, "N.clus:", N.clus, "clus.minNum:", clus.minNum))
    
    if(N.clus > clus.minNum){ 
      N.clus.prev <- N.clus
      clus.update.prev <- clus.update
    }else if(N.clus == clus.minNum){
      message(paste("i:", i, "ari.clus.resol:", ari.clus.resol.list[i], "N.clus:", N.clus, "clus.minNum:", clus.minNum ))
      break;
      
    }else if(N.clus < clus.minNum){
      if(i!=1){
        N.clus.next <- N.clus
        list.N.clus <- c(N.clus.prev, N.clus.next)
        w <- which.min(abs(list.N.clus - clus.minNum))
        if(w==1){clus.update <- clus.update.prev}
      }
      break;
    }
  }
  
  obj$seurat_clusters = clus.update
  Idents(obj) = clus.update
  return(obj)
}

generate_clustering_results.minSize.obj <- function(obj.name, subtype, fixed.resol=fixed.resol.value){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(subtype))
  #***** 4. Identification of highly variable features (feature selection)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose =FALSE)
  #***** 5. Scaling the data 
  # (,so that highly-expressed genes do not dominate)
  
  tryCatch({
    test <- ScaleData(obj, features = rownames(obj), verbose =FALSE)},
    error = function(e) print("ScaleData - initial trial returns error...Let's try again!"),
    warning = function(w) print("ScaleData - initial trial returns warning...Let's try again!"),
    finally = NULL)
  
  obj <- ScaleData(obj, features = rownames(obj), verbose =FALSE)
  
  #***** 6. Perform linear dimensional reduction (PCA)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE)
  
  #***** 7. Cluster the cells
  N.PCs <- length(obj[['pca']])
  obj <- FindNeighbors(obj, dims = 1:N.PCs, k.param=20, n.trees=50, verbose =FALSE) ## (Jaccard similarity)
  
  ari.clus.resol.list = fixed.resol #c(1.2)
  N = length(ari.clus.resol.list)
  message(paste("ARI.clus.update:TRUE --- try different resolution values... in [1.2 -0.1, 0.01 ]"))
  
  for(i in 1:N){
    tmp <- FindClusters(obj, resolution = ari.clus.resol.list[i], verbose =FALSE) ##previously, resolution = 0.5, apply modularity optimization techniques to iteratively group cells together 
    clus.update = tmp$seurat_clusters
    N.clus <- length(table(clus.update))
    message(paste("i:",i, "N.clus:", N.clus, "clus.minNum:", clus.minNum))
    
    if(N.clus > clus.minNum){ 
      N.clus.prev <- N.clus
      clus.update.prev <- clus.update
    }else if(N.clus == clus.minNum){
      message(paste("i:", i, "ari.clus.resol:", ari.clus.resol.list[i], "N.clus:", N.clus, "clus.minNum:", clus.minNum ))
      break;
      
    }else if(N.clus < clus.minNum){
      if(i!=1){
        N.clus.next <- N.clus
        list.N.clus <- c(N.clus.prev, N.clus.next)
        w <- which.min(abs(list.N.clus - clus.minNum))
        if(w==1){clus.update <- clus.update.prev}
      }
      break;
    }
  }
  
  
  obj$seurat_clusters = clus.update
  Idents(obj) = clus.update
  clusSize <- table(obj$seurat_clusters)
  clus.minSize = 0.01
  min.clusSize <- round(sum(clusSize)*clus.minSize) #1000
  w.smallClus <- which(clusSize <min.clusSize)
  N.smallClus <- length(w.smallClus)
  
  if( N.smallClus > 0 ){
    message(paste("apply min.clusSize:", min.clusSize))
    Idents(obj) <- ARI.clus_update.snn(obj, w.smallClus)
    obj$seurat_clusters <- Idents(obj)
    message("updated clustering results! --- assinged the minor cells to the nearest neighbor cell groups.")
  }
  
  # use snn matrix to get the same N.groups as real subgroups
  ARI.clus_update.snn <- function(obj, w.smallClus){
    # cluster level to remove
    smallClus <- names(table(Idents(obj)))[w.smallClus]
    
    ## update the cluster info
    #:  Assign the minor cell groups (small size clusters) to the nearest neighbor cells.
    message("generate neighbor cells clus info matrix...")
    
    #neighbor.info = obj@neighbors$RNA.nn@nn.idx
    snn.matrix = obj@graphs$RNA_snn
    snn.list <-data.frame(snn.matrix)#lapply(seq_len(ncol(snn.matrix)), function(i) snn.matrix[,i])
    clus.org.list <- as.character(Idents(obj))
    clus.remove.list <- lapply(clus.org.list,function(i){unique(c(i, smallClus))})
    
    clus.list <- rep(list(clus.org.list), length(snn.list))
    snames.list <- rep(list(names(Idents(obj))), length(snn.list))
    
    up.fun<- function(snn, clus.org, clus, snames){
      snn<-unlist(snn)
      snn.order = order(snn, decreasing=TRUE)
      snn.order.clus <- clus[snn.order]
      #w.next <- match(FALSE, grepl(clus.org, snn.order.clus))
      clus.org.txt <-paste("[",paste(clus.org, collapse="|"),"]")
      w.next <- match(FALSE, grepl(clus.org.txt, snn.order.clus))
      clus.update<- snn.order.clus[w.next]
      return(clus.update)
    }
    
    #clus.update=mapply(up.fun, snn=snn.list, clus.org=clus.org.list, clus=clus.list, snames=snames.list)
    clus.update=mapply(up.fun, snn=snn.list, clus.org=clus.remove.list, clus=clus.list, snames=snames.list)
    names(clus.update) = names(Idents(obj))
    
    clusToUpdate = data.frame(clus.org=as.character(Idents(obj)), clus.update=clus.update)
    
    N.smallClus <- length(w.smallClus)
    for(i in 1: N.smallClus){
      w<-which(clusToUpdate$clus.org==smallClus[i])
      Idents(obj)[w] <- clusToUpdate$clus.update[w]
    }
    
    return(Idents(obj))
  }
  
  return(obj)
}


# updated to have ARI, silhouette, NMI info in the object
generate_clustering_results.Vresol.obj <- function(obj.name, subtype){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(subtype))
  #***** 4. Identification of highly variable features (feature selection)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose =FALSE)
  #***** 5. Scaling the data 
  # (,so that highly-expressed genes do not dominate)
  
  tryCatch({
    test <- ScaleData(obj, features = rownames(obj), verbose =FALSE)},
    error = function(e) print("ScaleData - initial trial returns error...Let's try again!"),
    warning = function(w) print("ScaleData - initial trial returns warning...Let's try again!"),
    finally = NULL)
  
  obj <- ScaleData(obj, features = rownames(obj), verbose =FALSE)
  
  #***** 6. Perform linear dimensional reduction (PCA)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE)
  
  #***** 7. Cluster the cells
  N.PCs <- length(obj[['pca']])
  if(N.PCs > 20){ N.PCs <- 20}
  obj <- FindNeighbors(obj, dims = 1:N.PCs, k.param=20, n.trees=50, verbose =FALSE) ## (Jaccard similarity)
  
  ari.clus.resol.list = c(seq(1.2,0.1, -0.2), 0.1, 0.01)
  N = length(ari.clus.resol.list)
  message(paste("ARI.clus.update:TRUE --- try different resolution values... in [1.2 -0.1, 0.01 ]"))
  
  for(i in 1:N){
    tmp <- FindClusters(obj, resolution = ari.clus.resol.list[i], verbose =FALSE) ##previously, resolution = 0.5, apply modularity optimization techniques to iteratively group cells together 
    clus.update = tmp$seurat_clusters
    N.clus <- length(table(clus.update))
    message(paste("i:",i, "N.clus:", N.clus, "clus.minNum:", clus.minNum))
    
    if(N.clus > clus.minNum){ 
      N.clus.prev <- N.clus
      clus.update.prev <- clus.update
    }else if(N.clus == clus.minNum){
      message(paste("i:", i, "ari.clus.resol:", ari.clus.resol.list[i], "N.clus:", N.clus, "clus.minNum:", clus.minNum ))
      break;
      
    }else if(N.clus < clus.minNum){
      if(i!=1){
        N.clus.next <- N.clus
        list.N.clus <- c(N.clus.prev, N.clus.next)
        w <- which.min(abs(list.N.clus - clus.minNum))
        if(w==1){clus.update <- clus.update.prev}
      }
      break;
    }
  }
  
  
  obj$seurat_clusters = clus.update
  Idents(obj) = clus.update

  obj@misc$ARI <- aricode::ARI(obj$seurat_clusters, subtype)
  obj@misc$NMI <- aricode::NMI(obj$seurat_clusters, subtype)
  
  #Compute silhouette score for each cluster
  # : The silhouette width varies between -1 and 1, 
  #   for which a negative score means that observations may probably be in the wrong cluster, 
  #   while positive score means that observations are well clustered.
  
  #Compute distance matrix to UMAP coordinates
  N.PCs <- length(obj[['pca']])
  if(N.PCs > 20){ N.PCs <- 20}
  obj <- RunUMAP(obj, dims = 1:N.PCs)
  distance_matrix <- dist(Embeddings(obj[['umap']])[, 1:2])
  clusters <- as.character(subtype)  
  c.list = names(table(clusters))
  N = length( names(table(clusters)))
  for(i in 1:N){ clusters[clusters==c.list[i]] <- i }
  sil <- cluster::silhouette(as.numeric(clusters), dist = distance_matrix)
  ss <-   summary(sil)
  obj@misc$ASW <- round(ss$avg.width, digits=5)

  
  return(obj)
}





# feature selection method HVG 
generate_clustering_results.Vresol.obj_features <- function(obj.name, subtype, method, nfeatures, dist=NULL){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(subtype))
  method = method
  dist = dist
  nfeatures = nfeatures
  
  #***** 4. Identification of highly variable features (feature selection)
  features <-featureSelection(obj.name = obj, method=method, nfeatures = nfeatures, dist= dist)
  Seurat::VariableFeatures(object = obj) <- features

  #***** 5. Scaling the data 
  # (,so that highly-expressed genes do not dominate)
  tryCatch({
    test <- ScaleData(obj, features = NULL, verbose =FALSE)},
    error = function(e) print("ScaleData - initial trial returns error...Let's try again!"),
    warning = function(w) print("ScaleData - initial trial returns warning...Let's try again!"),
    finally = NULL)
  
  # done for each feature 
  obj <- ScaleData(obj, features = NULL, verbose =FALSE)
  
  
  #***** 6. Perform linear dimensional reduction (PCA)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE)
  
  #***** 7. Cluster the cells
  N.PCs <- length(obj[['pca']])
  if(N.PCs > 20){ N.PCs <- 20}
  obj <- FindNeighbors(obj, dims = 1:N.PCs, k.param=20, n.trees=50, verbose =FALSE) ## (Jaccard similarity)

  ari.clus.resol.list = c(seq(1.2,0.1, -0.2), 0.1, 0.01)
  N = length(ari.clus.resol.list)
  message(paste("ARI.clus.update:TRUE --- try different resolution values... in [1.2 -0.1, 0.01 ]"))

  for(i in 1:N){
    tmp <- FindClusters(obj, resolution = ari.clus.resol.list[i], verbose =FALSE) ##previously, resolution = 0.5, apply modularity optimization techniques to iteratively group cells together
    clus.update = tmp$seurat_clusters
    N.clus <- length(table(clus.update))
    message(paste("i:",i, "N.clus:", N.clus, "clus.minNum:", clus.minNum))

    if(N.clus > clus.minNum){
      N.clus.prev <- N.clus
      clus.update.prev <- clus.update
    }else if(N.clus == clus.minNum){
      message(paste("i:", i, "ari.clus.resol:", ari.clus.resol.list[i], "N.clus:", N.clus, "clus.minNum:", clus.minNum ))
      break;

    }else if(N.clus < clus.minNum){
      if(i!=1){
        N.clus.next <- N.clus
        list.N.clus <- c(N.clus.prev, N.clus.next)
        w <- which.min(abs(list.N.clus - clus.minNum))
        if(w==1){clus.update <- clus.update.prev}
      }
      break;
    }
  }


  obj$seurat_clusters = clus.update
  Idents(obj) = clus.update

  obj@misc$ARI <- aricode::ARI(obj$seurat_clusters, subtype)
  obj@misc$NMI <- aricode::NMI(obj$seurat_clusters, subtype)
  
  #Compute silhouette score for each cluster
  # : The silhouette width varies between -1 and 1, 
  #   for which a negative score means that observations may probably be in the wrong cluster, 
  #   while positive score means that observations are well clustered.
  
  #Compute distance matrix to UMAP coordinates
  N.PCs <- length(obj[['pca']])
  if(N.PCs > 20){ N.PCs <- 20}
  obj <- RunUMAP(obj, dims = 1:N.PCs)
  distance_matrix <- dist(Embeddings(obj[['umap']])[, 1:2])
  clusters <- as.character(subtype)  
  c.list = names(table(clusters))
  N = length( names(table(clusters)))
  for(i in 1:N){ clusters[clusters==c.list[i]] <- i }
  sil <- cluster::silhouette(as.numeric(clusters), dist = distance_matrix)
  ss <-   summary(sil)
  obj@misc$ASW <- round(ss$avg.width, digits=5)

  return(obj)
}

# feature selection method HVG 
featureSelection <- function(obj.name, method, nfeatures = 2000, dist=NULL){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  method = method
  dist = dist
  nfeatures = nfeatures
  
  #***** 4. Identification of highly variable features (feature selection)
  if(method == "HVG.vst"){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures, verbose =FALSE)
  }
  
  if(method == "HVG.mvp"){
    obj <- FindVariableFeatures(obj, selection.method = "mvp", nfeatures = nfeatures, verbose =FALSE)
  }
  
  if(method == "HVG.disp"){
    obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = nfeatures, verbose =FALSE)
  }
  
  if(method == "HDG"){
    #https://cran.r-project.org/web/packages/glmpca/index.html
    #The objective function to be minimized is simply the overall deviance
    #library(glmpca)
    library(scry)
    countData<- obj[["RNA"]]@counts
    
    if(dist=="poisson"){
      deviances <- scry::devianceFeatureSelection(countData, fam=c("poisson"))
    }else if(dist== "binomial"){
      deviances <- scry::devianceFeatureSelection(countData, fam=c("binomial")) # default farm is binomial
    }
    
    o <- order(deviances, decreasing=TRUE)
    HDG <- rownames(countData)[o][1:nfeatures]
    features <- HDG
    
    Seurat::VariableFeatures(object = obj) <- HDG
  } 
  
  
  if(method == "HEG.cpm"){
    x<- obj[["RNA"]]@counts
    cpm <- edgeR::cpm(x)
    genemeans<- rowMeans(cpm)
    o <- order(genemeans, decreasing=TRUE)
    HEG <- rownames(cpm)[o][1:nfeatures]
    Seurat::VariableFeatures(object = obj)  <- HEG
  } 
  
  if(method == "HVG.plus.HDG"){
    N.allF = nrow(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N.allF, verbose =FALSE)
    temp.hvg <-  Seurat::VariableFeatures(object = obj)
    
    library(scry)
    countData<- obj[["RNA"]]@counts
    deviances <- scry::devianceFeatureSelection(countData, fam=c("poisson"))
    o <- order(deviances, decreasing=TRUE)
    temp.hdg <- rownames(countData)[o][1:N.allF]
    
    
    # need to check how many common genes are
    for (i in 200:N.allF){
      common.features <- intersect(temp.hvg[1:i], temp.hdg[1:i])
      N.cf <- length(common.features)
      if(N.cf >= nfeatures ){ 
        common.features = common.features[1:nfeatures]
        message(paste("i:",i,"#common.features:", length(common.features)))
        break; 
      }
    }
    print(paste("#number.hvg.plus.hdg:", length(common.features)))
    Seurat::VariableFeatures(object = obj) <- common.features
  }
  
  
  if(method == "HVG.plus.HEG"){
    N.allF = nrow(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N.allF, verbose =FALSE)
    temp.hvg <-  Seurat::VariableFeatures(object = obj)
    
    x<- obj[["RNA"]]@counts
    cpm <- edgeR::cpm(x)
    genemeans<- rowMeans(cpm)
    o <- order(genemeans, decreasing=TRUE)
    temp.heg <- rownames(cpm)[o][1:N.allF]
    
    # need to check how many common genes are
    for (i in 200:N.allF){
      common.features <- intersect(temp.hvg[1:i], temp.heg[1:i])
      N.cf <- length(common.features)
      if(N.cf >= nfeatures ){ 
        common.features = common.features[1:nfeatures]
        message(paste("i:",i,"#common.features:", length(common.features)))
        break; 
      }
    }
    print(paste("#number.hvg.plus.heg:", length(common.features)))
    Seurat::VariableFeatures(object = obj) <- common.features
  }
  
  
  if(method == "NBDrop"){
    # M3Drop - Michaelis-Menten Modelling of Dropouts for scRNASeq
    # dropout-based feature selection method
    
    library(M3Drop)
    counts<- obj[["RNA"]]@counts
    #norm<- obj[["RNA"]]@data
    
    # M3Drop requires a normalized but not log-transformed expression matrix,
    norm <- M3DropConvertData(counts, is.counts=TRUE)
    
    # M3Drop requires a normalized but not log-transformed expression matrix,
    # thus the above function can optionally de-log a log2-normalized expression matrix from any normalization method.
    #norm <- M3DropConvertData(obj[["RNA"]]@data, is.log=TRUE, pseudocount=1)
    
    # spikes <- sample(1:nrow(counts), 50);
    # spikes <- rownames(norm)[spikes];
    # spikes <- rownames(norm) %in% spikes;
    Features_consensus <- Consensus_FS(counts, norm, include_cors = FALSE); #is.spike=spikes, 
    #head(Features_consensus)
    geneNames = rownames(Features_consensus)
    
    NBDrop = geneNames[order(Features_consensus$DANB_drop)][1:nfeatures]
    NBDisp = geneNames[order(Features_consensus$DANB_var)][1:nfeatures]
    M3Drop = geneNames[order(Features_consensus$M3Drop)][1:nfeatures]
    
    Seurat::VariableFeatures(object = obj) <- NBDrop
  }
  
  if(method == "NBDisp"){
    # M3Drop - Michaelis-Menten Modelling of Dropouts for scRNASeq
    # dropout-based feature selection method
    
    library(M3Drop)
    counts<- obj[["RNA"]]@counts
    #norm<- obj[["RNA"]]@data
    
    # M3Drop requires a normalized but not log-transformed expression matrix,
    norm <- M3DropConvertData(counts, is.counts=TRUE)
    
    # M3Drop requires a normalized but not log-transformed expression matrix,
    # thus the above function can optionally de-log a log2-normalized expression matrix from any normalization method.
    #norm <- M3DropConvertData(obj[["RNA"]]@data, is.log=TRUE, pseudocount=1)
    
    # spikes <- sample(1:nrow(counts), 50);
    # spikes <- rownames(norm)[spikes];
    # spikes <- rownames(norm) %in% spikes;
    Features_consensus <- Consensus_FS(counts, norm,  include_cors = FALSE); #is.spike=spikes,
    #head(Features_consensus)
    geneNames = rownames(Features_consensus)
    
    NBDrop = geneNames[order(Features_consensus$DANB_drop)][1:nfeatures]
    NBDisp = geneNames[order(Features_consensus$DANB_var)][1:nfeatures]
    M3Drop = geneNames[order(Features_consensus$M3Drop)][1:nfeatures]
    
    Seurat::VariableFeatures(object = obj) <- NBDisp
  }
  
  if(method == "M3Drop"){
    # M3Drop - Michaelis-Menten Modelling of Dropouts for scRNASeq
    # dropout-based feature selection method
    
    library(M3Drop)
    counts<- obj[["RNA"]]@counts
    #norm<- obj[["RNA"]]@data
    
    # M3Drop requires a normalized but not log-transformed expression matrix,
    norm <- M3DropConvertData(counts, is.counts=TRUE)
    
    # M3Drop requires a normalized but not log-transformed expression matrix,
    # thus the above function can optionally de-log a log2-normalized expression matrix from any normalization method.
    #norm <- M3DropConvertData(obj[["RNA"]]@data, is.log=TRUE, pseudocount=1)
    
    # spikes <- sample(1:nrow(counts), 50);
    # spikes <- rownames(norm)[spikes];
    # spikes <- rownames(norm) %in% spikes;
    Features_consensus <- Consensus_FS(counts, norm,  include_cors = FALSE); #is.spike=spikes,
    #head(Features_consensus)
    geneNames = rownames(Features_consensus)
    
    NBDrop = geneNames[order(Features_consensus$DANB_drop)][1:nfeatures]
    NBDisp = geneNames[order(Features_consensus$DANB_var)][1:nfeatures]
    M3Drop = geneNames[order(Features_consensus$M3Drop)][1:nfeatures]
    
    Seurat::VariableFeatures(object = obj) <- M3Drop
  }
  
  
  if(method == "HIPPO"){
    library(HIPPO)
    countData<- obj[["RNA"]]@counts
    
    #X = SingleCellExperiment::counts(sce)
    subX = as.matrix(countData) #X
    subdf = HIPPO:::preprocess_heterogeneous(subX)
    df <- HIPPO:::compute_test_statistic(subdf)
    # zero proportion value
    #features = subdf[subdf$zvalue > z_threshold, ]
    o <- order(df$zvalue, decreasing=TRUE)
    hippo.F <- rownames(df)[o][1:nfeatures]
    features <- hippo.F
    
    Seurat::VariableFeatures(object = obj) <- hippo.F
  } 
  
  if(method == "DUBStepR"){
    #1.logN normalized data (using all genes) =>DUBStepR: 546 features 
    #  - optimise.features = T ->only 546 genes? (obj: 19464 features across 20000 samples)
    #  - optimise.features = F ->only 568 genes? (obj: 19464 features across 20000 samples)
    #2.scT.v2 normalized data(using filterd genes by scT.v2) =>DUBStepR: 68 features
    #3.logN normalized data(using filtered genes by scT.v2) =>DUBStepR: 16 features
    
    # library(DUBStepR)
    #dubstepR.out.opt <- DUBStepR::DUBStepR(input.data = obj@assays$RNA@data, min.cells = 0, optimise.features = F, k = 10, num.pcs = 20, error = 0)
    
    # dubstepR.out.opt <- DUBStepR::DUBStepR(input.data = obj@assays$RNA@data, min.cells = 0.05*ncol(obj), optimise.features = T, k = 10, num.pcs = 20, error = 0)
    # length(dubstepR.out.opt$optimal.feature.genes)
    # DUBStepR.opt.F <- dubstepR.out.opt$optimal.feature.genes[1:nfeatures]

    library(DUBStepR.TEST)
    dubstepR.out <- DUBStepR.TEST::DUBStepR.nFeatures(nFeatures=nfeatures, input.data = obj@assays$RNA@data,  min.cells = 0.05*ncol(obj), optimise.features = F, k = 10, num.pcs = 20, error = 0)
    DUBStepR.F = dubstepR.out$corr.info$feature.genes
    print(length(DUBStepR.F))
    
    features <- DUBStepR.F
    Seurat::VariableFeatures(object = obj) <- DUBStepR.F
  } 
  
  features <-  Seurat::VariableFeatures(object = obj)
  
  return(features)
}



generate_clustering_results.Vresol.obj_features_DUBStepR <- function( model.name, dataName,obj.name, subtype, method, nfeatures, dist=NULL, filterN="filter0"){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(subtype))
  method = method
  dist = dist
  nfeatures = nfeatures
  
  
  #***** 4. Identification of highly variable features (feature selection)
  #features <-featureSelection(obj.name = obj, method=method, nfeatures = nfeatures, dist= dist)
  
  dataName = dataName
  model.name = model.name
  if(filterN=="filter0"){file.dir="~/8.data/allFeatures_DUBStepR_filter0/"}
  if(filterN=="filter0.05"){file.dir="~/8.data/allFeatures_DUBStepR_filter0.05/"}
  features_filename = paste(file.dir,"Features.ordered.DUBStepR.", dataName,".", model.name, ".nF",nfeatures,".txt", sep="")
  DUBStepR.nF.ordered = read.table(features_filename)
  DUBStepR.nF.ordered = DUBStepR.nF.ordered[,1]
  
  features = DUBStepR.nF.ordered
  Seurat::VariableFeatures(object = obj) <- features
  
  message(paste("nF:", length(DUBStepR.nF.ordered), "-",   length(features)))
  
  #***** 5. Scaling the data 
  # (,so that highly-expressed genes do not dominate)
  tryCatch({
    test <- ScaleData(obj, features = NULL, verbose =FALSE)},
    error = function(e) print("ScaleData - initial trial returns error...Let's try again!"),
    warning = function(w) print("ScaleData - initial trial returns warning...Let's try again!"),
    finally = NULL)
  
  # done for each feature 
  obj <- ScaleData(obj, features = NULL, verbose =FALSE)
  
  
  #***** 6. Perform linear dimensional reduction (PCA)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE)
  
  #***** 7. Cluster the cells
  N.PCs <- length(obj[['pca']])
  if(N.PCs > 20){ N.PCs <- 20}
  obj <- FindNeighbors(obj, dims = 1:N.PCs, k.param=20, n.trees=50, verbose =FALSE) ## (Jaccard similarity)

  ari.clus.resol.list = c(seq(1.2,0.1, -0.2), 0.1, 0.01)
  N = length(ari.clus.resol.list)
  message(paste("ARI.clus.update:TRUE --- try different resolution values... in [1.2 -0.1, 0.01 ]"))

  for(i in 1:N){
    tmp <- FindClusters(obj, resolution = ari.clus.resol.list[i], verbose =FALSE) ##previously, resolution = 0.5, apply modularity optimization techniques to iteratively group cells together
    clus.update = tmp$seurat_clusters
    N.clus <- length(table(clus.update))
    message(paste("i:",i, "N.clus:", N.clus, "clus.minNum:", clus.minNum))

    if(N.clus > clus.minNum){
      N.clus.prev <- N.clus
      clus.update.prev <- clus.update
    }else if(N.clus == clus.minNum){
      message(paste("i:", i, "ari.clus.resol:", ari.clus.resol.list[i], "N.clus:", N.clus, "clus.minNum:", clus.minNum ))
      break;

    }else if(N.clus < clus.minNum){
      if(i!=1){
        N.clus.next <- N.clus
        list.N.clus <- c(N.clus.prev, N.clus.next)
        w <- which.min(abs(list.N.clus - clus.minNum))
        if(w==1){clus.update <- clus.update.prev}
      }
      break;
    }
  }


  obj$seurat_clusters = clus.update
  Idents(obj) = clus.update

  obj@misc$ARI <- aricode::ARI(obj$seurat_clusters, subtype)
  obj@misc$NMI <- aricode::NMI(obj$seurat_clusters, subtype)
  
  
  #Compute silhouette score for each cluster
  # : The silhouette width varies between -1 and 1, 
  #   for which a negative score means that observations may probably be in the wrong cluster, 
  #   while positive score means that observations are well clustered.
  
  #Compute distance matrix to UMAP coordinates
  N.PCs <- length(obj[['pca']])
  if(N.PCs > 20){ N.PCs <- 20}
  obj <- RunUMAP(obj, dims = 1:N.PCs)
  distance_matrix <- dist(Embeddings(obj[['umap']])[, 1:2])
  clusters <- as.character(subtype)  
  c.list = names(table(clusters))
  N = length( names(table(clusters)))
  for(i in 1:N){ clusters[clusters==c.list[i]] <- i }
  sil <- cluster::silhouette(as.numeric(clusters), dist = distance_matrix)
  ss <-   summary(sil)
  obj@misc$ASW <- round(ss$avg.width, digits=5)
  
  return(obj)
}


featuresExpGroupsV2 <- function( expGroups , features = VariableFeatures(obj)){
  EGR = list()
  EGR$G1 = length(which(features %in% expGroups$G1))/length(features)
  EGR$G2 = length(which(features %in% expGroups$G2))/length(features)
  EGR$G3 = length(which(features %in% expGroups$G3))/length(features)
  EGR$G4 = length(which(features %in% expGroups$G4))/length(features)
  EGR$G5 = length(which(features %in% expGroups$G5))/length(features)
  
  #EGR = sprintf("%.2f", EGR)
  return(EGR)

}


trueDEGRV2 <- function(degData, features = VariableFeatures(obj) ){
  DEGR = list()
  DEGR$all = length(which(features %in% degData$all))/(length(features))
  DEGR$up = length(which(features %in% degData$up_genes))/(length(features))
  DEGR$down = length(which(features %in% degData$down_genes))/(length(features))
  return(DEGR)
}


trueDEGRV3 <- function(degData, features = VariableFeatures(obj) ){
  DEGR = list()
  
  DEGR$all = length(which(features %in% degData$all))/(length(features))
  DEGR$up = length(which(features %in% degData$up_genes))/(length(features))
  DEGR$down = length(which(features %in% degData$down_genes))/(length(features))
  
  TPR
  FPR
  AUR #auROC()
  
  return(DEGR)
}
