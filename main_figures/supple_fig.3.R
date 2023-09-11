#supple. fig.3 
# : (5 subtypes clusters) lung_myeloid_scd_v3 in scTv2 results
# Top4(HDG, HEG, HIPPO, NBDrop)
# Low4(HVG.disp, HVG.vst, NBDisp, M3Drop)

#source("~/3.scripts/main_figures/supple_fig.3.R")

#generate_supple_fig.3(dataName="lung_myeloid_scd_v3", filter_rate=0)


generate_supple_fig.3 <- function(dataName, filter_rate){
  ####### data loading ################################
  # dir.result = "tmp.sim.DrHiDepth_featureNumbers"
  # file.path.result = paste("/store/juok/research/Research__scNormalization/results/", dir.result,  "/",sep = "")
  # runName="DrHiDepth_featureNumbers"
  # 
  ############# load real data ######################
  #n = as.numeric(DataN)
  filter_rate = filter_rate
  
  #DrHiDepthSimData v3&v4
  path.data = "/store/juok/research/Research__scNormalization/data/Dr.Hi.simData/"
  dataName = "lung_myeloid_scd_v3"
  
  file.list = dir(path.data)
  w = grep(paste(dataName,".data.",sep=""), file.list)
  datafile.list = file.list[w]

  library("SingleCellExperiment")
  fileName = paste(path.data, datafile.list, sep="")
  simData <- readRDS(fileName) 
  count.data = counts(simData)
  real.subtype <- simData$Group   
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
  count.data <- count.data[w, ]
  cpm.data.org = cpm.data[w,]
  
  # create seurat obj
  library(Seurat)
  N.minCells = 0
  obj.org <- Seurat::CreateSeuratObject(counts = count.data, project = "real", min.cells = N.minCells)
  obj.org[["percent.mt"]] <- PercentageFeatureSet(obj.org, pattern = "^MT-")
  message("created seurat obj...")
  obj.org
  
  method.list.top4 = c("scTv2.HDG", "scTv2.HEG.cpm",  "scTv2.HIPPO", "scTv2.NBDrop")
   
  pTop.list = generate_integrated_EGR.violin_umap_top4(obj=obj.org, real.subtype, cpm.data.org, method.list = method.list.top4, width=1)
  
  p.list1 = pTop.list$p.EGR.list
  p.list2 = pTop.list$p.umap.list
  x= length(p.list1)
  y=length(p.list2)
  for(i in 1:x){
    p1Name = paste("pT", i,sep="")
    assign(p1Name, p.list1[[i]])
  }
  
  for(j in 1:y){
    p2Name = paste("pTu", j, sep="")
    assign(p2Name, p.list2[[j]])
  }
  
  method.list.bottom4 = c("scTv2.HVG.disp", "scTv2.HVG.vst",  "scTv2.NBDisp", "scTv2.M3Drop")
  
  pBottom.list = generate_integrated_EGR.violin_umap_bottom4(obj=obj.org, real.subtype, cpm.data.org, method.list=method.list.bottom4,width=0.1)
  
  p.list1 = pBottom.list$p.EGR.list
  p.list2 = pBottom.list$p.umap.list
  x= length(p.list1)
  y=length(p.list2)
  for(i in 1:x){
    p1Name = paste("pB", i,sep="")
    assign(p1Name, p.list1[[i]])
  }
  
  for(j in 1:y){
    p2Name = paste("pBu", j, sep="")
    assign(p2Name, p.list2[[j]])
  }
  
  # for ARI & deg plot of all methods
  p<- generate_ARI.plot_sfig3(dataName.list = c("lung_myeloid_scd_v3"), models.list = c("scTv2"))

  p.ari <- p$p.ari
  p.deg <- p$p.deg
  
  # for ASW plot of all methods
  p.asw <- generate_ASW.plot(dataName.list= c("lung_myeloid_scd_v3"), models.list=c("scTv2"))
  
  # for N2S plot 
  p.n2s <- generate_N2S.plot(obj=obj.org, cpm.data.org)
  
  
  
  #******* summary.result plot
  # integrated EGR plot
  file.path.summary = "/store/juok/research/Research__scNormalization/results/summary.plot/"
  timestamp = paste0("___",format(Sys.time(), "%Y%m%d_%H%M%S"))
  plotName = paste("Supple_Fig.3__Top4_and_Bottom4_",dataName, "_", filter_rate, "_",timestamp,".pdf", sep="")
  
  pdfName = paste(file.path.summary, plotName, sep="")
  
  default.txt = paste("Supple Fig.3: Top 4 and Bottom 4 feature selections - ",  dataName, sep="")
                          # " - HighlyExp", 100*(1-filter_rate),"%", sep="")
  
  pdf(pdfName, width=23, height=25)
  library(gridExtra)
  library(grid)
  title.txt = default.txt
  
  gl = list(pT1,pT2,pT3,pT4,
            pTu1, pTu2,pTu3,pTu4,
            pB1,pB2,pB3,pB4,
            pBu1,pBu2, pBu3, pBu4, 
            p.ari, p.asw, p.deg, p.n2s)
  
  grid.arrange(
    grobs = gl,
    widths = c(1, 1, 1, 1),
    layout_matrix = rbind(c(1, 2, 3, 4),
                          c(5, 6, 7, 8),
                          c(9,10,11,12),
                          c(13,14,15,16),
                          c(17, 18, 19, 20)),
    top = textGrob(title.txt, gp = gpar(fontsize = 25, fontface="bold")) 
  )
  
  dev.off()
  
  message(paste("generated summary plot!..\n",pdfName))
  
  # p.umap.legend<- generate_UMAP.realSubtype_obj_features(model.name = method.list.bottom4[j],obj.name = obj.org, real.subtype = real.subtype, method="HDG", nfeatures = 500, dist="poisson", LPos="bottom")
  # pdf("legend_umap_lung_myeloid_scd_v3_fig.3.pdf", 20,10)
  # p.umap.legend
  # dev.off()
  
}



generate_integrated_EGR.violin_umap_top4 <-function(obj, real.subtype, cpm.data.org, method.list, width){
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
  
  p.EGR.list = list()
  p.umap.list = list()
  
  #2. scTv2.HDG
  obj.scTv2.HDG = obj.scTv2.HVG.vst
  i=1; j=1
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HDG, method="HDG", dist="poisson",width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j], obj.name = obj.scTv2.HDG, real.subtype = real.subtype, method="HDG", nfeatures = 500, dist="poisson")
  p.umap.list[[j]] = p
  print("method top1.scT v2.HDG...ARI Done!")
  
  #3.scTv2.HEG.cpm
  obj.scTv2.HEG.cpm = obj.scTv2.HVG.vst
  i=1; j=2
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HEG.cpm, method="HEG.cpm", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.HEG.cpm, real.subtype = real.subtype, method="HEG.cpm", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method top2.scT v1.HEG.cpm...ARI Done!")
  
  #17. "scTv2.HIPPO"
  obj.scTv2.HIPPO = obj.scTv2.HVG.vst
  i=1; j=3
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HIPPO, method="HIPPO", dist=NULL,width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.HIPPO, real.subtype = real.subtype, method="HIPPO", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method top 3.scTv2.HIPPO...ARI Done!")
  
  #6. "scTv2.NBDrop"
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
  p.list$p.EGR.list = p.EGR.list
  p.list$p.umap.list = p.umap.list
  
  return(p.list)
}

generate_integrated_EGR.violin_umap_bottom4 <-function(obj, real.subtype, cpm.data.org, method.list, width){
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
  
  p.EGR.list = list()
  p.umap.list = list()
  
  #21. scT v2.HVG.disp
  obj.scTv2.HVG.disp = obj.scTv2.HVG.vst
  i=1; j=1
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HVG.disp, method="HVG.disp", dist=NULL, width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j], obj.name = obj.scTv2.HVG.disp, real.subtype = real.subtype, method="HVG.disp", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method bottom 1.scT v2.HVG disp...ARI Done!")
  
  #3.scTv2.HVG.vst
  i=1; j=2
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.HVG.vst, method="HVG.vst", dist=NULL, width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.HVG.vst, real.subtype = real.subtype, method="HVG.vst", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method bottom 2.scT v2.HVG vst...ARI Done!")
  
  #17. "scTv2.NBDisp"
  obj.scTv2.NBDisp = obj.scTv2.HVG.vst
  i=1; j=3
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.NBDisp, method="NBDisp", dist=NULL, width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.NBDisp, real.subtype = real.subtype, method="NBDisp", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method bottom 3.scTv2.NBDisp..ARI Done!")
  
  #6. "scTv2.M3Drop"
  obj.scTv2.M3Drop = obj.scTv2.HVG.vst
  i=1; j=4
  p <- generate_EGR.violin.plot(model.name= method.list[j], cpm.data.org, obj.name= obj.scTv2.M3Drop, method="M3Drop", dist=NULL, width=width)
  p.EGR.list[[j]] = p
  p<- generate_UMAP.realSubtype_obj_features(model.name = method.list[j],obj.name = obj.scTv2.M3Drop, real.subtype = real.subtype, method="M3Drop", nfeatures = 500, dist=NULL)
  p.umap.list[[j]] = p
  print("method bottom 4.scTv2.M3Drop...ARI Done!")
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
    ggtitle(paste(model.name, "\n", "ASW:",ASW)) +
    theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
          axis.title = element_text(size=15, face="bold"),
          plot.title = element_text(size=15, hjust=0.5, face="bold"),
          legend.position=LPos) 
          #legend.position="none"  
          # legend.title = element_text(size=14),
          # legend.text=element_text(size=12))

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

  #legend.title = element_text(size=14),
  #legend.text=element_text(size=12))
  
  #+ geom_jitter(shape=16, position=position_jitter(0.2)) + geom_boxplot(width=0.1)
  
  #p
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
  N.ylim = 30#range(N2S.gmean)[2]
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


generate_ASW.plot <- function(dataName.list = c("lung_myeloid_scd_v3"), models.list = c("scTv2")){
  library("ggplot2")
  library(RColorBrewer)
  library(scales)
  
  models.list = models.list
  dir.result = "tmp.sim.DrHiDepth_featureNumbers"
  file.path.result = paste("/store/juok/research/Research__scNormalization/results/", dir.result,  "/",sep = "")
  runName="DrHiDepth_featureNumbers"
  dataName.list = c("lung_myeloid_scd_v3")
  
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
      
      tableName = paste("DrHiDepth",dataName.list[n], filter_rate, sep="_")
      ASWfileName = paste(file.path.result, "ASW.table", tableName,".txt", sep="")
      
      models.list = c("scTv2") #c( "scTv2","logN")
      N.m = length(models.list)
      # k for model
      for(k in 1:N.m){
        # s for plot style-(ARI, EGR, DEG)
        s=1
        # ARI plot
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


generate_ARI.plot_sfig3 <- function(dataName.list = c("lung_myeloid_scd_v3"), models.list = c("scTv2")){
  library("ggplot2")
  library(RColorBrewer)
  library(scales)
  
  models.list = models.list
  dir.result = "tmp.sim.DrHiDepth_featureNumbers"
  file.path.result = paste("/store/juok/research/Research__scNormalization/results/", dir.result,  "/",sep = "")
  runName="DrHiDepth_featureNumbers"
  
  dataName.list = c("lung_myeloid_scd_v3")
  
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
      
      tableName = paste("DrHiDepth",dataName.list[n], filter_rate, sep="_")
      ARIfileName = paste(file.path.result, "ARI.table", tableName,".txt", sep="")
      DEGallfileName = paste(file.path.result, "DEGall.table", tableName,".txt", sep="")
      
      models.list = c("scTv2") #c( "scTv2","logN")
      N.m = length(models.list)
      # k for model
      for(k in 1:N.m){
        # s for plot style-(ARI, EGR, DEG)
       for(s in 1:2){
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
          
          pName = paste("p",s,k, sep="")
          print(paste("p.s.k:",pName,sep=""))
          assign(pName, p)
        }#s1
        
        # DEG rate - all plot
        if(s==2){ 
          dfileName.list = c(DEGallfileName)
          dytitle.list = c("All DEG True Rate")
          N.deg= 1
          for(d in 1:N.deg){
            dfileName = dfileName.list[d]
            dytitle = dytitle.list[d]
            if(file.exists(dfileName)){
              tmp = read.table(dfileName, quote="", sep="\t") 
              colnames(tmp) =  c("scTv2.HVG.vst","scTv2.HDG", "scTv2.HEG.cpm",
                                 "scTv2.HVG.plus.HDG", "scTv2.HVG.plus.HEG",
                                 "scTv2.NBDrop", "scTv2.NBDisp", "scTv2.M3Drop", 
                                 "scTv2.HIPPO","scTv2.HVG.disp","scTv2.DUBStepR",
                                 "logN.HVG.vst", "logN.HDG", "logN.HEG.cpm",
                                 "logN.HVG.plus.HDG", "logN.HVG.plus.HEG",
                                 "logN.NBDrop", "logN.NBDisp", "logN.M3Drop",
                                 "logN.HIPPO",  "logN.HVG.disp", "logN.DUBStepR")
            }else{ 
              message("d:",d, dfileName, "  result file do not exist!" )
            }
            
            m = colnames(tmp)
            #models.list = c(scTv2","logN")
            w <- grep(paste(models.list[k],sep=""), m)
            print(paste("k:", k, " - model:", models.list[k]))
            score.table = tmp[,w]
            
            # score plots
            plot.table = as.matrix(score.table)
            method.list = rep(gsub(paste(models.list[k],".",sep=""), "",colnames(plot.table)), each=nrow(plot.table))
            nfeatures.list = rep(rownames(plot.table), ncol(plot.table))
            score.list = c(plot.table) #sprintf("%.3f",c(plot.table))
            plot.data = data.frame(method=method.list, nfeatures=nfeatures.list, trueDEGR=score.list)
            plot.data$nfeatures = factor(plot.data$nfeatures, levels=rownames(plot.table))
            l=c("HDG", "HEG.cpm", "HIPPO", "DUBStepR", "HVG.vst", "HVG.plus.HEG", "HVG.plus.HDG", "HVG.disp","NBDrop", "NBDisp", "M3Drop")
            #f=paste(models.list[k], l, sep=".")
            plot.data$method = factor(plot.data$method, levels=l)
            ptitle = models.list[k]
            
            p <- ggplot(plot.data,  aes(x=nfeatures, y=trueDEGR, color=method, group=method, linetype=method, size=method)) + 
              #geom_point(col=3) + 
              geom_line() + coord_cartesian(ylim = c(0.47, 0.60)) +
              ggtitle(ptitle) +
              ylab(dytitle) + theme_bw() + 
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
            
            pName = paste("p", s, k, d, sep="")
            print(paste("p.s.k.d:",pName,sep=""))
            assign(pName, p)
          }#d for deg all up down
        }#s2 for style, ari , deg
        
       }#s  
      }#k for model 
    }#f for filter
  }#n for dataN
  
  p.list = list()
  p.list$p.ari = p11
  p.list$p.deg = p211
  return(p.list)
}
