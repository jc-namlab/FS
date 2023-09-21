# supple. fig.4 -> (facs10) logN version
# : top4(ARI,EGR umap * scTv2, logN) -> top4(3*2 plot) + bottom4 (3*2 plot)in 1 plot 

#source("~/3.scripts/main_figures/supple_fig.4.R")

#generate_supple_fig.4(dataName.list="Zheng_sorted", filter_rate=0)

# Top4(HDG, HEG, HIPPO, NBDrop)
# Low4(HVG.disp, HVG.vst, NBDisp, M3Drop)

generate_supple_fig.4 <- function(dataName.list, filter_rate){
  ####### data loading ################################
  # dir.result = "tmp.sim.DrHiDepth_featureNumbers"
  # file.path.result = paste("/store/juok/research/Research__scNormalization/results/", dir.result,  "/",sep = "")
  # runName="DrHiDepth_featureNumbers"
  # 
  
  library("DuoClustering2018")
  library("SingleCellExperiment")
  
  # n = as.numeric(DataN)
  # filter_rate = filter_rate
  source("/store/juok/research/Research__scNormalization/scripts/generate_RealData_Abdelaal_Benchmark_datasets.R")
  
  ############# load real data ######################
  #n = as.numeric(DataN)
  filter_rate = 0#filter_rate
  dataName.list = c("Zheng_sorted") #facs10
  
  data = generate_RealData_Abdelaal_Benchmark(n=1)
  count.data = data$countData
  real.subtype = data$realSubtype
  rownames(count.data)=gsub("_","-", rownames(count.data))
  
  # all zero genes filtering
  rowsum = rowSums(count.data)
  w.rm = which(rowsum==0)
  if(length(w.rm)>0){ count.data = count.data[-w.rm, ]}
  
  # cpm exp level 
  cpm.data = edgeR::cpm(count.data)
  gmean = rowMeans(cpm.data)
  cutoff = quantile(gmean, as.numeric(filter_rate))
  w <- which(gmean >= cutoff)
  count.data <- count.data[w,]
  cpm.data.org = cpm.data[w,]
  
  
  # create seurat obj
  library(Seurat)
  N.minCells = 0
  obj.org <- CreateSeuratObject(counts = count.data, project = "real", min.cells = N.minCells)
  obj.org[["percent.mt"]] <- PercentageFeatureSet(obj.org, pattern = "^MT-")
  message("created seurat obj...")
  obj.org
  
  # for ARI plot of all methods
  p.ari<- generate_ARI.plot_fig3(dataName.list = c("Zheng_sorted"), models.list = c("logN"))
  
  # for ASW plot of all methods
  p.asw <- generate_ASW.plot(dataName.list= c("Zheng_sorted"), models.list=c("logN"))
  
  # for N2S plot 
  p.n2s <- generate_N2S.plot(obj=obj.org, cpm.data.org)
  
  # for top 4 panel
  models.list = c("scTv2", "logN")
  N.m= length(models.list)
  for(m in 2:N.m){
    # if(m == 1){
    #   method.list.top4 = c("logN.HDG", "logN.HEG.cpm",  "logN.HIPPO", "logN.NBDrop")
    #   pTop.list = generate_integrated_ARI_EGR.violin_umap_scTv2_top4(obj=obj.org, real.subtype, cpm.data.org, method.list=method.list.top4, width=1)
    # }
    if(m == 2){
      method.list.top4 = c("logN.HDG", "logN.HEG.cpm",  "logN.HIPPO", "logN.NBDrop")
      pTop.list = generate_integrated_ARI_EGR.violin_umap_logN_top4(obj=obj.org, real.subtype, cpm.data.org, method.list=method.list.top4, width=1)
    }

    p.list2 = pTop.list$p.EGR.list
    p.list3 = pTop.list$p.umap.list
    
    #egr
    y=length(p.list2)
    for(j in 1:y){
      p2Name = paste("pT2",m, j, sep="")
      assign(p2Name, p.list2[[j]])
    }
    
    #umap
    z=length(p.list3)
    for(k in 1:z){
      p3Name = paste("pT3",m, k, sep="")
      assign(p3Name, p.list3[[k]])
    }
  }# m for models

  # for bottom 4 panel
  models.list = c("scTv2", "logN")
  N.m= length(models.list)
  for(m in 2:N.m){
    # if(m == 1){
    #   method.list.bottom4 = c("scTv2.HVG.disp", "scTv2.HVG.vst",  "scTv2.NBDisp", "scTv2.M3Drop")
    #   pBottom.list = generate_integrated_ARI_EGR.violin_umap_scTv2_bottom4(obj=obj.org, real.subtype, cpm.data.org, method.list=method.list.bottom4, width=0.1)
    # }
    if(m == 2){
      method.list.bottom4 = c("logN.HVG.disp", "logN.HVG.vst",  "logN.NBDisp", "logN.M3Drop")
      pBottom.list = generate_integrated_ARI_EGR.violin_umap_logN_bottom4(obj=obj.org, real.subtype, cpm.data.org, method.list=method.list.bottom4, width=0.1)
    }

    p.list2 = pBottom.list$p.EGR.list
    p.list3 = pBottom.list$p.umap.list

    #egr
    y=length(p.list2)
    for(j in 1:y){
      p2Name = paste("pB2",m, j, sep="")
      assign(p2Name, p.list2[[j]])
    }
    
    #umap
    z=length(p.list3)
    for(k in 1:z){
      p3Name = paste("pB3",m, k, sep="")
      assign(p3Name, p.list3[[k]])
    }
  }# m for models
  
  gl.all = list(p.ari,  p.asw, p.n2s,              
              pT221, pT222,pT223,pT224, 
              pT321, pT322,pT323,pT324,
              pB221, pB222,pB223,pB224,
              pB321, pB322,pB323,pB324)
  
  generate_all_pdf(gl=gl.all, dataName="facs10", filter_rate, w=20, h=20)
  
  
  # p.umap.legend<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.logN.M3Drop, real.subtype = real.subtype, method="M3Drop", nfeatures = 500, dist=NULL, LPos="bottom")
  # pdf("legend_umap_facs10_fig.3.pdf", 20,10)
  # p.umap.legend
  # dev.off()

}

generate_all_pdf <- function(gl, dataName, filter_rate, w=20, h=20){
  #******* summary.result plot
  # integrated EGR plot
  file.path.summary = "/store/juok/research/Research__scNormalization/results/summary.plot/"
  timestamp = paste0("___",format(Sys.time(), "%Y%m%d_%H%M%S"))
  plotName = paste("Supple_fig.4__Top4_and_Bottom4_",dataName, "_", filter_rate, "_",timestamp,".pdf", sep="")
  
  pdfName = paste(file.path.summary, plotName, sep="")
  
  default.txt = paste("Supple Fig.4: Top4 and Bottom4 ", "feature selections in ",  dataName, sep="")
                                        # " - HighlyExp", 100*(1-filter_rate),"%", sep="")
  
  pdf(pdfName, width=w, height=h)
  library(gridExtra)
  library(grid)
  title.txt = default.txt
  
  gl = list(p.ari,  p.asw, p.n2s,              
            pT221, pT222,pT223,pT224, 
            pT321, pT322,pT323,pT324,
            pB221, pB222,pB223,pB224,
            pB321, pB322,pB323,pB324)

  grid.arrange(
    nrow= 5,
    grobs = gl,
    widths = c(1,1,1,1), heights = c(1.5, 1,1,1,1),
    layout_matrix = rbind(c(1, 2, NA, 3),
                          c(4, 5, 6, 7),
                          c(8, 9, 10, 11),
                          c(12,13, 14, 15),
                          c(16, 17,18, 19)),
    top = textGrob(title.txt, gp = gpar(fontsize = 25, fontface="bold"))
  )

  dev.off()
  
  message(paste("generated all pdf!..\n",pdfName))
  
} 


generate_integrated_ARI_EGR.violin_umap_scTv2_top4 <-function(obj, real.subtype, cpm.data.org, method.list, width){
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(sctransform)
  library(Linnorm)
  library(ggplot2)
  library(aricode)
  library(edgeR)
  library(DESeq2)
  library(scran)
  library(cluster)
  library("glmGamPoi")
  library("sctransformV2Test")
  library("ggplot2")
  
  #***** Run different normalization method
  
  #1.scT.v2 HVG.counts
  obj.scTv2.HVG.vst = obj
  scT.v2_out <- sctransformV2Test::vst_S2N(obj[["RNA"]]@counts, vst.flavor = "v2")
  obj.scTv2.HVG.vst[["RNA"]]@data <- scT.v2_out$y
  
  genesToUse = rownames(scT.v2_out$y)
  cellsToUse = colnames(scT.v2_out$y)
  obj[["RNA"]]@counts = obj[["RNA"]]@counts[genesToUse, cellsToUse]
  obj[["RNA"]]@data = obj[["RNA"]]@data[genesToUse, cellsToUse]
  obj.scTv2.HVG.vst[["RNA"]]@counts = obj[["RNA"]]@counts

  nfeatures.list = c(500, 1000, 2000, 4000) 
  N.f = length(nfeatures.list)
  ARI.mat = matrix(0, ncol=4, nrow=4)
  colnames(ARI.mat) =  method.list
  rownames(ARI.mat) = nfeatures.list
  
  p.EGR.list = list()
  p.umap.list = list()

  #1. scTv2.HDG
  obj.scTv2.HDG = obj.scTv2.HVG.vst
  i=1; j=1
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HDG, method="HDG", dist="poisson",width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j], obj.name = obj.scTv2.HDG, real.subtype = real.subtype, method="HDG", nfeatures = 500, dist="poisson")
  p.umap.list[[j]] = p
  print("method top1.scTv2.HDG...ARI Done!")
  
  #2.scTv2.HEG.cpm
  obj.scTv2.HEG.cpm = obj.scTv2.HVG.vst
  i=1; j=2
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HEG.cpm, method="HEG.cpm", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.HEG.cpm, real.subtype = real.subtype, method="HEG.cpm", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method top2.scTv2.HEG.cpm...ARI Done!")
  
  #3. "scTv2.HIPPO"
  obj.scTv2.HIPPO = obj.scTv2.HVG.vst
  i=1; j=3
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HIPPO, method="HIPPO", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.HIPPO, real.subtype = real.subtype, method="HIPPO", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method top 3.scTv2.HIPPO...ARI Done!")
  
  #4. "scTv2.NBDrop"
  obj.scTv2.NBDrop = obj.scTv2.HVG.vst
  i=1; j=4
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.NBDrop, method="NBDrop", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.NBDrop, real.subtype = real.subtype, method="NBDrop", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method top 4.scTv2.NBDrop...ARI Done!")
  #######################################
  
  message("finish! return a p.list ")
  
  p.list = list()
  p.list$ARI.mat = ARI.mat
  p.list$p.EGR.list = p.EGR.list
  p.list$p.umap.list = p.umap.list
  
  return(p.list)
}

generate_integrated_ARI_EGR.violin_umap_logN_top4 <-function(obj, real.subtype, cpm.data.org, method.list, width){
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(sctransform)
  library(Linnorm)
  library(ggplot2)
  library("DuoClustering2018")
  library("SingleCellExperiment")
  library(aricode)
  library(edgeR)
  library(DESeq2)
  library(scran)
  library(cluster)
  library("glmGamPoi")
  library("sctransformV2Test")
  library("ggplot2")
  
  
  #***** Run different normalization method
  
  #1.scT.v2 HVG.counts
  obj.scTv2.HVG.vst = obj
  scT.v2_out <- sctransformV2Test::vst_S2N(obj[["RNA"]]@counts, vst.flavor = "v2")
  obj.scTv2.HVG.vst[["RNA"]]@data <- scT.v2_out$y
  
  genesToUse = rownames(scT.v2_out$y)
  cellsToUse = colnames(scT.v2_out$y)
  obj[["RNA"]]@counts = obj[["RNA"]]@counts[genesToUse, cellsToUse]
  obj[["RNA"]]@data = obj[["RNA"]]@data[genesToUse, cellsToUse]
  obj.scTv2.HVG.vst[["RNA"]]@counts = obj[["RNA"]]@counts
  
  obj.logN.HVG.vst = obj
  obj.logN.HVG.vst <- NormalizeData(obj.logN.HVG.vst, normalization.method = "LogNormalize", scale.factor = 10000)
  
  
  nfeatures.list = c(500, 1000, 2000, 4000) 
  N.f = length(nfeatures.list)

  p.EGR.list = list()
  p.umap.list = list()
  
  #1. scTv2.HDG
  obj.logN.HDG = obj.logN.HVG.vst
  i=1; j=1
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.HDG, method="HDG", dist="poisson",width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j], obj.name = obj.logN.HDG, real.subtype = real.subtype, method="HDG", nfeatures = 500, dist="poisson", LPos="none")
  p.umap.list[[j]] = p
  print("method top1.logN.HDG...ARI Done!")
  
  #2.logN.HEG.cpm
  obj.logN.HEG.cpm = obj.logN.HVG.vst
  i=1; j=2
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.HEG.cpm, method="HEG.cpm", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.logN.HEG.cpm, real.subtype = real.subtype, method="HEG.cpm", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method top2.logN.HEG.cpm...ARI Done!")
  
  #3. "logN.HIPPO"
  obj.logN.HIPPO = obj.logN.HVG.vst
  i=1; j=3
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.HIPPO, method="HIPPO", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.logN.HIPPO, real.subtype = real.subtype, method="HIPPO", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method top 3.logN.HIPPO...ARI Done!")
  
  #4. "logN.NBDrop"
  obj.logN.NBDrop = obj.logN.HVG.vst
  i=1; j=4
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.NBDrop, method="NBDrop", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.logN.NBDrop, real.subtype = real.subtype, method="NBDrop", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method top 4.logN.NBDrop...ARI Done!")
  #######################################
  
  message("finish! return a p.list ")
  
  p.list = list()
  p.list$p.EGR.list = p.EGR.list
  p.list$p.umap.list = p.umap.list
  
  return(p.list)
}

generate_integrated_ARI_EGR.violin_umap_scTv2_bottom4 <-function(obj, real.subtype, cpm.data.org, method.list, width){
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(sctransform)
  library(Linnorm)
  library(ggplot2)
  library("DuoClustering2018")
  library("SingleCellExperiment")
  library(aricode)
  library(edgeR)
  library(DESeq2)
  library(scran)
  library(cluster)
  library("glmGamPoi")
  library("sctransformV2Test")
  library("ggplot2")
  
  
  #***** Run different normalization method
  
  #1.scT.v2 HVG.counts
  obj.scTv2.HVG.vst = obj
  scT.v2_out <- sctransformV2Test::vst_S2N(obj[["RNA"]]@counts, vst.flavor = "v2")
  obj.scTv2.HVG.vst[["RNA"]]@data <- scT.v2_out$y
  
  genesToUse = rownames(scT.v2_out$y)
  cellsToUse = colnames(scT.v2_out$y)
  obj[["RNA"]]@counts = obj[["RNA"]]@counts[genesToUse, cellsToUse]
  obj[["RNA"]]@data = obj[["RNA"]]@data[genesToUse, cellsToUse]
  obj.scTv2.HVG.vst[["RNA"]]@counts = obj[["RNA"]]@counts
  
  nfeatures.list = c(500, 1000, 2000, 4000) 
  N.f = length(nfeatures.list)
  ARI.mat = matrix(0, ncol=4, nrow=4)
  colnames(ARI.mat) =  method.list
  rownames(ARI.mat) = nfeatures.list
  
  p.EGR.list = list()
  p.umap.list = list()
  
  #1. scT v2.HVG.disp
  obj.scTv2.HVG.disp = obj.scTv2.HVG.vst
  i=1; j=1
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HVG.disp, method="HVG.disp", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j], obj.name = obj.scTv2.HVG.disp, real.subtype = real.subtype, method="HVG.disp", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 1.scT v2.HVG disp...ARI Done!")
  
  #2.scTv2.HVG.vst
  i=1; j=2
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HVG.vst, method="HVG.vst", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.HVG.vst, real.subtype = real.subtype, method="HVG.vst", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 2.scT v2.HVG vst...ARI Done!")
  
  #3. "scTv2.NBDisp"
  obj.scTv2.NBDisp = obj.scTv2.HVG.vst
  i=1; j=3
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.NBDisp, method="NBDisp", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.NBDisp, real.subtype = real.subtype, method="NBDisp", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 3.scTv2.NBDisp..ARI Done!")
  
  #4. "scTv2.M3Drop"
  obj.scTv2.M3Drop = obj.scTv2.HVG.vst
  i=1; j=4
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.M3Drop, method="M3Drop", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.M3Drop, real.subtype = real.subtype, method="M3Drop", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 4.scTv2.M3Drop...ARI Done!")
  #######################################
  
  message("finish! return a p.list ")
  
  p.list = list()
  p.list$ARI.mat = ARI.mat
  p.list$p.EGR.list = p.EGR.list
  p.list$p.umap.list = p.umap.list
  
  return(p.list)
}

generate_integrated_ARI_EGR.violin_umap_logN_bottom4 <-function(obj, real.subtype, cpm.data.org, method.list, width){
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(sctransform)
  library(Linnorm)
  library(ggplot2)
  library("DuoClustering2018")
  library("SingleCellExperiment")
  library(aricode)
  library(edgeR)
  library(DESeq2)
  library(scran)
  library(cluster)
  library("glmGamPoi")
  library("sctransformV2Test")
  library("ggplot2")
  
  
  #***** Run different normalization method
  
  #1.scT.v2 HVG.counts
  obj.scTv2.HVG.vst = obj
  scT.v2_out <- sctransformV2Test::vst_S2N(obj[["RNA"]]@counts, vst.flavor = "v2")
  obj.scTv2.HVG.vst[["RNA"]]@data <- scT.v2_out$y
  
  genesToUse = rownames(scT.v2_out$y)
  cellsToUse = colnames(scT.v2_out$y)
  obj[["RNA"]]@counts = obj[["RNA"]]@counts[genesToUse, cellsToUse]
  obj[["RNA"]]@data = obj[["RNA"]]@data[genesToUse, cellsToUse]
  obj.scTv2.HVG.vst[["RNA"]]@counts = obj[["RNA"]]@counts
  
  obj.logN.HVG.vst = obj
  obj.logN.HVG.vst <- NormalizeData(obj.logN.HVG.vst, normalization.method = "LogNormalize", scale.factor = 10000)
  
  nfeatures.list = c(500, 1000, 2000, 4000) 
  N.f = length(nfeatures.list)
 
  p.EGR.list = list()
  p.umap.list = list()
  
  #1. scT v2.HVG.disp
  obj.logN.HVG.disp = obj.logN.HVG.vst
  i=1; j=1
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.HVG.disp, method="HVG.disp", dist=NULL ,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j], obj.name = obj.logN.HVG.disp, real.subtype = real.subtype, method="HVG.disp", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 1.logN.HVG disp...ARI Done!")
  
  #2.logN.HVG.vst
  i=1; j=2
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.HVG.vst, method="HVG.vst", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.logN.HVG.vst, real.subtype = real.subtype, method="HVG.vst", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 2.logN HVG vst...ARI Done!")
  
  #3. "logN.NBDisp"
  obj.logN.NBDisp = obj.logN.HVG.vst
  i=1; j=3
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.NBDisp, method="NBDisp", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.logN.NBDisp, real.subtype = real.subtype, method="NBDisp", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 3.logN.NBDisp..ARI Done!")
  
  #4. "logN.M3Drop"
  obj.logN.M3Drop = obj.logN.HVG.vst
  i=1; j=4
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.logN.M3Drop, method="M3Drop", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.logN.M3Drop, real.subtype = real.subtype, method="M3Drop", nfeatures = 500, dist=NULL, LPos="none")
  p.umap.list[[j]] = p
  print("method bottom 4.logN.M3Drop...ARI Done!")
  #######################################
  
  message("finish! return a p.list ")
  
  p.list = list()
  p.list$p.EGR.list = p.EGR.list
  p.list$p.umap.list = p.umap.list
  
  return(p.list)
}

# feature selection method HVG 
generate_UMAP.realSubtype_obj_features <- function(model.name, obj.name, real.subtype, method, nfeatures, dist=NULL, LPos="none"){
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(real.subtype))
  method = method
  dist = dist
  nfeatures = nfeatures
  
  # add meta data - cell info, cell attr, batch info..
  names(real.subtype) <- colnames(x = obj)
  
  obj <- AddMetaData(object = obj, metadata = real.subtype, col.name='real_subtypes')
  #head(obj@meta.data)
  
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
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose =FALSE, seed.use=42)
  
  #Compute distance matrix to UMAP coordinates
  N.PCs <- length(obj[['pca']])
  if(N.PCs > 20){ N.PCs <- 20}
  obj <- RunUMAP(obj, dims = 1:N.PCs)

  
  # ASW
  #Compute distance matrix to UMAP coordinates
  distance_matrix <- dist(Embeddings(obj[['umap']])[, 1:2])
  
  clusters <- as.character(real.subtype)  
  c.list = names(table(clusters))
  N = length( names(table(clusters)))
  for(i in 1:N){ clusters[clusters==c.list[i]] <- i }
  sil <- cluster::silhouette(as.numeric(clusters), dist = distance_matrix)
  ss <-   summary(sil)
  ASW <- round(ss$avg.width, digits=5)
  ###################
  library("ggplot2")
  library("scales")
  col.cd4.v1 = c("mediumblue", "#6baed6",  "#2171b5", "#00B0F6"  )
  col.cd4.v2 = c( "#39B600", "#00B0F6",  "#2171b5", "mediumblue" )
  col.cd4.v3 = c("mediumblue", "#6baed6", "#39B600", "#00B0F6"  )
  col.cd4.v4 = c("mediumblue", "#6baed6", "#00BFC4", "#00B0F6"  )
  
  c1=scales::hue_pal()(10)
  mycols = c(c1[1:3], col.cd4.v4, c1[8:10])
  #pt.size=1
  p.umap <- DimPlot(obj, reduction ="umap", pt.size=0.1, group.by="real_subtypes", cols=mycols) + 
    #ggtitle(paste(model.name, "-",nfeatures, "features", "\n", "ASW:",ASW)) +
    ggtitle(paste(model.name, "\n", "ASW:",ASW)) +
    theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
          axis.title = element_text(size=15, face="bold"),
          plot.title = element_text(size=15, hjust=0.5, face="bold"),
          legend.position=LPos) 
          #legend.position="none"
          #legend.title = element_text(size=14),
          #legend.text=element_text(size=12))
          
          #theme(legend.position="bottom")
  
  return(p.umap)
  
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
    # dubstepR.out.opt <- DUBStepR::DUBStepR(input.data = obj@assays$RNA@data, min.cells = 0.05*ncol(obj), optimise.features = T, k = 10, num.pcs = 20, error = 0)
    # length(dubstepR.out.opt$optimal.feature.genes)
    # DUBStepR.opt.F <- dubstepR.out.opt$optimal.feature.genes[1:nfeatures]
    
    # library(DUBStepR.TEST)
    # dubstepR.out <- DUBStepR.TEST::DUBStepR.minNumGenes(minNumGenes = nfeatures, input.data = obj@assays$RNA@data,  min.cells = 0, optimise.features = F, k = 10, num.pcs = 20, error = 0)
    # DUBStepR.F = dubstepR.out$corr.info$feature.genes
    # print(length(DUBStepR.F))
    
    features <- DUBStepR.F
    Seurat::VariableFeatures(object = obj) <- DUBStepR.F
  } 
  
  features <-  Seurat::VariableFeatures(object = obj)
  
  return(features)
}

generate_EGR.violin.plot <- function(model.name, cpm.data.org, obj.name, method, dist, width=1){
  nfeatures.list = c(500, 1000, 2000, 4000) 
  nF = length(nfeatures.list)
  EGR.cpm.exp.levels = list() 
  EGR.nfeatures = list()
  
  for(f in 1:nF){
    # #***** 4. Identification of highly variable features (feature selection)
    features <-featureSelection(obj.name, method=method, nfeatures = nfeatures.list[f], dist=dist)
    Seurat::VariableFeatures(object = obj.name) <- features
    
    # cpm exp level 
    #cpm.data = edgeR::cpm(count.data)
    gmean.org = rowMeans(cpm.data.org)
    gmean.org.ordered = sort(gmean.org)
    #w=which(names(gmean.org) %in% features)
    #fgmean = gmean.org[w]
    #EGR.cpm.exp.levels[[f]] = fgmean
    fgmean.rank = which(names(gmean.org.ordered) %in% features)
    EGR.cpm.exp.levels[[f]] = fgmean.rank
    EGR.nfeatures[[f]] = rep(nfeatures.list[f], as.numeric(nfeatures.list[f]))
  }
  
  df = data.frame( nFeatures = unlist(EGR.nfeatures),
                   cpm.exp.levels = unlist(EGR.cpm.exp.levels))
  head(df)
  
  df$nFeatures = factor(df$nFeatures, levels=nfeatures.list)
  
  #plot - EGR.violin
  library("ggplot2")
  library(RColorBrewer)
  library(scales)
  # > range(gmean)
  # [1] 9.886796e-03 5.435183e+03
  #q = quantile(gmean, c(0.2,0.4,0.6,0.8))
  # 20%        40%        60%        80% 
  # 8.684053  29.447334  66.551586 144.925399
  
  N.ylim = nrow(cpm.data.org)
  p <- ggplot(df, aes(x=nFeatures, y=cpm.exp.levels, fill=nFeatures)) +
    scale_fill_grey() + 
    geom_violin(width=width) + coord_cartesian(ylim = c(0, N.ylim)) + 
    ggtitle(model.name) + theme_bw() + 
    labs( y = "cpm.exp.levels(Rank)") +
    theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
          axis.title = element_text(size=15, face="bold"),
          plot.title = element_text(size=15, hjust=0.5, face="bold"),
          legend.position="none")

  #+ geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1)
  #legend.title = element_text(size=14),
  #legend.text=element_text(size=12))
  #p
  return(p)
}

# feature selection method HVG 
generate_clustering_results.Vresol.obj_features <- function(obj.name, subtype, method, nfeatures, dist=NULL){
  library(aricode)
  library(Seurat)
  library(sctransform)
  #library(normTEST)
  #library(normT)
  library(Linnorm)
  library(edgeR)
  library(DESeq2)
  library(scran)
  library(cluster)
  #library(DeNorm)
  library("glmGamPoi")
  library("sctransformV2Test")
  library(M3Drop)
  library(HIPPO)
  library(DUBStepR.TEST)
  
  
  #We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. 
  # ari.clus.resol = 1.5 # >1.2 
  obj = obj.name
  clus.minNum=length(table(subtype))
  method = method
  dist = dist
  
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
  
  # #Compute silhouette score for each cluster
  # # : The silhouette width varies between -1 and 1, 
  # #   for which a negative score means that observations may probably be in the wrong cluster, 
  # #   while positive score means that observations are well clustered.
  # 
  # #Compute distance matrix to UMAP coordinates
  # N.PCs <- length(obj[['pca']])
  # if(N.PCs > 20){ N.PCs <- 20}
  # obj <- RunUMAP(obj, dims = 1:N.PCs)
  # distance_matrix <- dist(Embeddings(obj[['umap']])[, 1:2])
  # clusters <- as.character(subtype)  
  # c.list = names(table(clusters))
  # N = length( names(table(clusters)))
  # for(i in 1:N){ clusters[clusters==c.list[i]] <- i }
  # sil <- cluster::silhouette(as.numeric(clusters), dist = distance_matrix)
  # ss <-   summary(sil)
  # obj@misc$ASW <- round(ss$avg.width, digits=5)
  
  return(obj)
}


generate_ARI.plot_fig3 <- function(dataName.list = c("Zheng_sorted"), models.list = c("logN")){
  library("ggplot2")
  library(RColorBrewer)
  library(scales)
  
  models.list = models.list
  dir.result = "tmp.real.Abdelaal.Benchmark_featureNumbers"
  file.path.result = paste("/store/juok/research/Research__scNormalization/results/", dir.result,  "/",sep = "")
  runName="RealAbdelaalBenchmark_featureNumbers"
  
  dataName.list = c("Zheng_sorted")

  #1. data 
  N.data = length(dataName.list)
  
  #2. dataName 
  for(n in 1:N.data){
    n = 1
    print(paste("n:",n))
    
    #3. filter_rate
    filter.list = c(0)#c(0, 0.1, 0.2) #, 0.5, 0.7)
    N.f = length(filter.list)
    for(f in 1:N.f){
      print(paste("f:",f))
      filter_rate = filter.list[f]
      
      tableName = paste("RealDataAbdelaal",dataName.list[n], filter_rate, sep="_")
      ARIfileName = paste(file.path.result, "ARI.table", tableName,".txt", sep="")
      
      models.list = c("logN") #c( "scTv2","logN")
      N.m = length(models.list)
      # k for model
      for(k in 1:N.m){
        # s for plot style-(ARI, EGR, DEG)
        s=1
        # ARI plot
        if(s==1){
          if(file.exists(ARIfileName)){
            tmp = read.table(ARIfileName, quote="", sep="\t") 
            colnames(tmp) =  c("scTv2.HVG.vst","scTv2.HDG", "scTv2.HEG.cpm",
                               "scTv2.HVG.plus.HDG", "scTv2.HVG.plus.HEG",
                               "scTv2.NBDrop", "scTv2.NBDisp", "scTv2.M3Drop", 
                               "scTv2.HIPPO","scTv2.HVG.disp","scTv2.DUBStepR",
                               "logN.HVG.vst", "logN.HDG", "logN.HEG.cpm",
                               "logN.HVG.plus.HDG", "logN.HVG.plus.HEG",
                               "logN.NBDrop", "logN.NBDisp", "logN.M3Drop",
                               "logN.HIPPO",  "logN.HVG.disp", "logN.DUBStepR")
          }else{ 
            message("n:",n, "  result file do not exist!" )
          }
          m = colnames(tmp)
          w <- grep(paste(models.list[k],sep=""), m)
          print(paste("k:", k, " - model:", models.list[k]))
          score.table = tmp[,w]
          
          # score plots
          plot.table = as.matrix(score.table)
          method.list = rep(gsub(paste(models.list[k],".",sep=""), "",colnames(plot.table)), each=nrow(plot.table))
          nfeatures.list = rep(rownames(plot.table), ncol(plot.table))
          score.list = c(plot.table) #sprintf("%.3f",c(plot.table))
            
          plot.data = data.frame(method=method.list, nfeatures=nfeatures.list, ARI=score.list)
          plot.data$nfeatures = factor(plot.data$nfeatures, levels=rownames(plot.table))
          l=c("HDG", "HEG.cpm", "HIPPO", "DUBStepR", "HVG.vst", "HVG.plus.HEG", "HVG.plus.HDG", "HVG.disp","NBDrop", "NBDisp", "M3Drop")
          #f=paste(models.list[k], l, sep=".")
          plot.data$method = factor(plot.data$method, levels=l)
          ptitle = models.list[k]
          
          p <- ggplot(plot.data,  aes(x=nfeatures, y=ARI, color=method, group=method, linetype=method, size=method)) + 
            geom_point(col=3) + 
            geom_line() + 
            ggtitle(ptitle) + theme_bw() + 
            theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
                  axis.title = element_text(size=15, face="bold"),
                  plot.title = element_text(size=15, hjust=0.5, face="bold"),
                  #legend.title = element_text(size=14),
                  #legend.text=element_text(size=12))
                  legend.position="none")
          
          st=plot_style(m=unique(plot.data$method))
          
          p = p + scale_linetype_manual(values = st$plot_lty)+
            scale_colour_manual(values = st$plot_col)+
            scale_size_manual(values = st$plot_size)+
            guides(guide_legend(nrow = 3, byrow = TRUE))
          
          # pName = paste("p",s,k, sep="")
          # print(paste("p.s.k:",pName,sep=""))
          # assign(pName, p)
        }#s1
      }#k for model 
    }#f for filter
  }#n for dataN
  
  return(p)
}

plot_style <- function(m){
  #display.brewer.all()
  #m = unique(plot.data$method)
  user_palette = c()
  user_palette[grep("DUBStepR",m)]<-brewer.pal(9,"Set1")[3]
  user_palette[grep("HDG",m)]<- brewer.pal(9,"Set1")[1]
  user_palette[grep("HEG.cpm",m)]<- brewer.pal(9,"Set1")[8] #4 pupple
  user_palette[grep("HIPPO",m)]<- brewer.pal(9,"Set1")[6] #8 pink 
  
  user_palette[grep("HVG.disp",m)]<- brewer.pal(9,"Blues")[5]
  user_palette[grep("HVG.plus.HEG",m)]<- brewer.pal(9,"Blues")[7]
  user_palette[grep("HVG.plus.HDG",m)]<- brewer.pal(9,"Blues")[6]
  user_palette[grep("HVG.vst",m)]<- brewer.pal(9,"Blues")[8]
  
  user_palette[grep("NBDrop",m)]<- brewer.pal(9,"Oranges")[6]
  user_palette[grep("M3Drop",m)]<- brewer.pal(9,"Oranges")[5]
  user_palette[grep("NBDisp",m)]<- brewer.pal(9,"Oranges")[4]
  
  plot_col= user_palette
  plot_lty= rep('solid',length(m)) #solid", "longdash", "dashed", "dotted"
  plot_size=rep(1,length(m))
  
  names(plot_col) = m
  names(plot_lty)=names(plot_size)=names(plot_col)
  
  #model=models.list[k]
  plot_lty[['HVG.plus.HDG']]=plot_lty[['HVG.plus.HEG']]="4121"
  plot_lty[['HVG.disp']]=plot_lty[['M3Drop']]="11" #'11'
  
  plot_size[[ "HVG.vst"]]=plot_size[['NBDrop']]=1.3
  plot_size[['HVG.disp']]=plot_size[['M3Drop']]=0.7
  
  
  # model=models.list[k]
  # plot_lty[[paste(model,'HVG.plus.HDG',sep=".")]]=plot_lty[[paste(model,'HVG.plus.HEG',sep=".")]]="4121"
  # plot_lty[[paste(model,'HVG.disp',sep=".")]]=plot_lty[[paste(model,'M3Drop',sep=".")]]="11" #'11'
  # 
  # plot_size[[paste(model, "HVG.vst",sep=".")]]=plot_size[[paste(model,'NBDrop',sep=".")]]=1.3
  # plot_size[[paste(model,'HVG.disp',sep=".")]]=plot_size[[paste(model,'M3Drop',sep=".")]]=0.7
  
  ######################################################################################### 
  
  plot_style = list()
  plot_style$plot_size = plot_size
  plot_style$plot_lty = plot_lty
  plot_style$plot_col = plot_col
  return(plot_style)
}

generate_ASW.plot <- function(dataName.list = c("Zheng_sorted"), models.list = c("logN")){
  library("ggplot2")
  library(RColorBrewer)
  library(scales)
  
  models.list = models.list
  dir.result = "tmp.real.Abdelaal.Benchmark_featureNumbers"
  file.path.result = paste("/store/juok/research/Research__scNormalization/results/", dir.result,  "/",sep = "")
  runName="RealAbdelaalBenchmark_featureNumbers"
  
  dataName.list = c("Zheng_sorted")
  
  #1. data 
  N.data = length(dataName.list)
  
  #2. dataName 
  for(n in 1:N.data){
    n = 1
    print(paste("n:",n))
    
    #3. filter_rate
    filter.list = c(0)#c(0, 0.1, 0.2) #, 0.5, 0.7)
    N.f = length(filter.list)
    for(f in 1:N.f){
      print(paste("f:",f))
      filter_rate = filter.list[f]
      
      tableName = paste("RealDataAbdelaal",dataName.list[n], filter_rate, sep="_")
      ASWfileName = paste(file.path.result, "ASW.table", tableName,".txt", sep="")
      
      models.list = c("logN") #c( "scTv2","logN")
      N.m = length(models.list)
      # k for model
      for(k in 1:N.m){
        # s for plot style-(ASW, EGR, DEG)
        s=1
        # ASW plot
        if(s==1){
          if(file.exists(ASWfileName)){
            tmp = read.table(ASWfileName, quote="", sep="\t") 
            colnames(tmp) =  c("scTv2.HVG.vst","scTv2.HDG", "scTv2.HEG.cpm",
                               "scTv2.HVG.plus.HDG", "scTv2.HVG.plus.HEG",
                               "scTv2.NBDrop", "scTv2.NBDisp", "scTv2.M3Drop", 
                               "scTv2.HIPPO","scTv2.HVG.disp","scTv2.DUBStepR",
                               "logN.HVG.vst", "logN.HDG", "logN.HEG.cpm",
                               "logN.HVG.plus.HDG", "logN.HVG.plus.HEG",
                               "logN.NBDrop", "logN.NBDisp", "logN.M3Drop",
                               "logN.HIPPO",  "logN.HVG.disp", "logN.DUBStepR")
          }else{ 
            message("n:",n, "  result file do not exist!" )
          }
          m = colnames(tmp)
          w <- grep(paste(models.list[k],sep=""), m)
          print(paste("k:", k, " - model:", models.list[k]))
          score.table = tmp[,w]
          
          # score plots
          plot.table = as.matrix(score.table)
          method.list = rep(gsub(paste(models.list[k],".",sep=""), "",colnames(plot.table)), each=nrow(plot.table))
          nfeatures.list = rep(rownames(plot.table), ncol(plot.table))
          score.list = c(plot.table) #sprintf("%.3f",c(plot.table))
          
          plot.data = data.frame(method=method.list, nfeatures=nfeatures.list, ASW=score.list)
          plot.data$nfeatures = factor(plot.data$nfeatures, levels=rownames(plot.table))
          l=c("HDG", "HEG.cpm", "HIPPO", "DUBStepR", "HVG.vst", "HVG.plus.HEG", "HVG.plus.HDG", "HVG.disp","NBDrop", "NBDisp", "M3Drop")
          #f=paste(models.list[k], l, sep=".")
          plot.data$method = factor(plot.data$method, levels=l)
          ptitle = models.list[k]
          
          p <- ggplot(plot.data,  aes(x=nfeatures, y=ASW, color=method, group=method, linetype=method, size=method)) + 
            geom_point(col=3) + 
            geom_line() + 
            ggtitle(ptitle) + theme_bw() + 
            theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
                  axis.title = element_text(size=15, face="bold"),
                  plot.title = element_text(size=15, hjust=0.5, face="bold"),
                  legend.title = element_text(size=14),
                  legend.text=element_text(size=12))
          #legend.position="none")
          
          st=plot_style(m=unique(plot.data$method))
          
          p = p + scale_linetype_manual(values = st$plot_lty)+
            scale_colour_manual(values = st$plot_col)+
            scale_size_manual(values = st$plot_size)+
            guides(guide_legend(nrow = 3, byrow = TRUE))
          
          # pName = paste("p",s,k, sep="")
          # print(paste("p.s.k:",pName,sep=""))
          # assign(pName, p)
        }#s1
      }#k for model 
    }#f for filter
  }#n for dataN
  
  return(p)
}

generate_N2S.plot <- function(obj=obj.org, cpm.data.org){
  
  # cpm exp level 
  gmean.all = rowMeans(cpm.data.org)
  q = quantile(gmean.all, c(0.25,0.5 ,0.75))
  w1 = which(gmean.all < q[1])
  w2 = which(gmean.all <q[2] & gmean.all >=q[1])
  w3 = which(gmean.all <q[3] & gmean.all >=q[2])
  w4 = which(gmean.all >= q[3])
  expGroups = list()
  expGroups$G1 = names(gmean.all)[w1]
  expGroups$G2 = names(gmean.all)[w2]
  expGroups$G3 = names(gmean.all)[w3]
  expGroups$G4 = names(gmean.all)[w4]

  scT.v2_out <- sctransformV2Test::vst_S2N(obj[["RNA"]]@counts, vst.flavor = "v2")
  N2S.gmean <- rowMeans(1/scT.v2_out$s2n)
  N2SG <- N2S_groups(expGroups, N2S.gmean)
  
  exp.level = c(rep("exp.G1", length(N2SG$G1)),
                rep("exp.G2", length(N2SG$G2)),
                rep("exp.G3", length(N2SG$G3)),
                rep("exp.G4", length(N2SG$G4)))
  N2S = c(N2SG$G1, N2SG$G2, N2SG$G3, N2SG$G4)
  features = c(names(N2SG$G1), names(N2SG$G2), names(N2SG$G3), names(N2SG$G4))
  
  data.boxplot = data.frame(features = features,
                            N2S = N2S,
                            exp.level = exp.level)
  ###############################################
  library(ggplot2)
  data.boxplot$exp.level <- factor(data.boxplot$exp.level, levels = c( "exp.G1","exp.G2", "exp.G3" , "exp.G4"))
  
  ptitle = ""
  N.ylim = range(N2S.gmean)[2]
  p <- ggplot(data.boxplot, aes(x=exp.level, y=N2S, color=exp.level)) +  geom_boxplot() + 
    #labs(title=paste("CPM Data(all features)")) + 
    coord_cartesian(ylim = c(0, N.ylim)) + 
    ggtitle(ptitle) + theme_bw() + 
    theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
          axis.title = element_text(size=15, face="bold"),
          plot.title = element_text(size=15, hjust=0.5, face="bold"),
          #legend.title = element_text(size=14),
          #legend.text=element_text(size=12),
          legend.position = "none")
  
    # dataMedian <-  dplyr::summarise(group_by(data.boxplot, exp.level), MD = as.numeric(sprintf("%.3f", median(N2S))))
    # p + geom_text(data = dataMedian, aes(exp.level, MD, label = MD),
    #         position = position_dodge(width = 0.8), size = 3, vjust = -0.5)

    # theme(axis.text =element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
    #       axis.title = element_text(size=15, face="bold"),
    #       plot.title = element_text(size=15, hjust=0.5, face="bold"),
    #       #legend.title = element_text(size=14),
    #       #legend.text=element_text(size=12),
    #       legend.position = "none")  + theme_bw()
  

  
  return(p)
  
}



N2S_groups <- function( expGroups , N2S.gmean){
  N2SG = list()
  N2SG$G1 = N2S.gmean[which(names(N2S.gmean) %in% expGroups$G1)]
  N2SG$G2 = N2S.gmean[which(names(N2S.gmean) %in% expGroups$G2)]
  N2SG$G3 = N2S.gmean[which(names(N2S.gmean) %in% expGroups$G3)]
  N2SG$G4 = N2S.gmean[which(names(N2S.gmean) %in% expGroups$G4)]
  
  return(N2SG)
}

