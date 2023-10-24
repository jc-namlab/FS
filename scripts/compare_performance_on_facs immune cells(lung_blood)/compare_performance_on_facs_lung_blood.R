{
  library(SummarizedExperiment)
  library(Seurat)
  library(magrittr)
  library(tidyverse)
  library(SeuratObject)
}
{
  #bulkSeq codes are originated from https://github.com/krasnowlab/HLCA/
  bulkSeq <- read.csv('weissman_blood_bulk_FPKM.csv', sep = "\t")
  rownames(bulkSeq) <- bulkSeq[,'gene_id']
  bulkSeq <- bulkSeq[,grep('FPKM',colnames(bulkSeq))]
  bulkSeq <- bulkSeq[-grep('^MIR', rownames(bulkSeq)),]
  bulkSeq <- bulkSeq[-grep('^SNO', rownames(bulkSeq)),]
  bulkSeq <- bulkSeq[-c(7969, 7970, 7971, 7972, 7973, 7974, 7976, 7977, 7978, 7979, 7980, 7982),]
  bulkSeq <- log1p(bulkSeq)
  
  # Reorder bulkSeq table to match SS2 order
  bulkSeq <- bulkSeq[,c(20,11,10,5,8,9,19,4,7,18,12,16,15,17,3,21,14,13,2,1,6)]
  
  # Remove bulk celltypes with low numbers of differentially expressed genes from other sequenced populations
  bulkSeq <- bulkSeq[,-c(6,13,20)]
  
  
  colData <- DataFrame(label.origin=colnames(bulkSeq),
                       label.refine=colnames(bulkSeq),
                       row.names=colnames(bulkSeq))
  
  # Combining naive, unswitched, and switched B for DE
  colData$label.refine[c(1,2,3)]=colData$label.refine[1]
  # Combining CD8 central memory and effector memory T for DE
  colData$label.refine[c(4,5)]=colData$label.refine[4]
  # Combining CD4 central memory and effector memory T for DE
  colData$label.refine[c(7,8)]=colData$label.refine[7]
  #   # Combining CD4 central memory and effector memory T for DE
  colData$label.refine[c(10,11)]=colData$label.refine[10]
  # colData$label.refine[c(2,3,5,8,11)]='none'
  
  
  ref.mrna<-SummarizedExperiment(assays=list(logcounts=bulkSeq),colData = colData)
  ref.mrna.use<-ref.mrna[,!str_detect(pattern='Eosinophil',string=colnames(ref.mrna))]
  ref.mrna.use$label.refine%<>%gsub(pattern='_FPKM',replacement='')
  data.use='SS2'
  for(patient in c("P1","P2","P3")){
    if(data.use=='SS2'){
      load(paste0('facs_normal_lung_blood_seurat_ntiss.',patient,'.anno.gencode.20200616.RC4.Robj'))
      assign(paste0("ntiss.",patient,".anno.gencode"),UpdateSeuratObject(get(paste0("ntiss.",patient,".anno.gencode"))))
    }else if(data.use=='10Xdroplet'){
      if(patient=="P2"){
        load(paste0('droplet_normal_lung_seurat_ntiss10x.',patient,'.anno.20191002.RC4.Robj'))
      }else{
        load(paste0('droplet_normal_lung_blood_seurat_ntiss10x.',patient,'.anno.20191002.RC4.Robj'))
      }
      assign(paste0("ntiss10x.",patient,".anno"),UpdateSeuratObject(get(paste0("ntiss10x.",patient,".anno"))))
    }
  }
}



immune_clus.use<-data.frame(orig.clus=c("B", "CD8+ Memory/Effector T", "CD8+ Naive T",
                                   "CD4+ Memory/Effector T", "CD4+ Naive T", "Natural Killer",
                                   "Neutrophil", "Basophil/Mast 1", "Plasmacytoid Dendritic",
                                   "Myeloid Dendritic Type 2", "Classical Monocyte", "Nonclassical Monocyte" ),
                            bulkseq_clus=c("Naïve_B_Cell","Central_Memory_CD8_T_Cell", "Naïve_CD8_T_Cell", "Central_Memory_CD4_T_Cell", "Naïve_CD4_T_Cell", "Immature_NK_Cell", "Neutrophil",
                                        "Basophil", "Plasmacytoid_Dendritic_Cell", "Myeloid_Dendritic_Cell", "Classical_Monocyte", "Non.Classical_Monocyte"))


{
  largecluster=T
  data.use='SS2'
  NPC=20
  dim.reduction='umap'
  
  
  `%ni%`=Negate(`%in%`)
  library(gplots)
  library(magrittr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
  library(BiocParallel)
  library(data.table)
  multicoreParam <- MulticoreParam(workers = 6)
  library(Seurat)
  library(celldex)
  library(SingleR)
  library(scales)
  library(gridExtra)
  library(grid)
  
  N_immune_gene='All'
  
  dir.base='/custom_output_dir/'
  dir.create(dir.base)
  
  dim.reduction='umap'
  
  dark_vivid_colors <- c( '#cc0000','#000054',"#0000FF", '#186b1b', "#66a999", '#ff6000',  '#c0927e', '#bfbfbf', '#531f82', '#c20b73', '#666300', '#FFBB34')#'darkblue'
  brighter_colors <- c('#FF9896', '#899CFF','#82EEFD','#4be176', '#8DD3C7', '#ffad66',  '#c0927e', '#bfbfbf','#8c6ab0', '#ffb2e8', '#d2e061','#ffe79e')
  
  dark_vivid_colors<-dark_vivid_colors[c(1:8,11)]
  brighter_colors<-brighter_colors[c(1:8,11)]
  
  color_extend<-rbind(dark_vivid_colors,brighter_colors)%>%as.vector()
}
if(data.use=='SS2'){
  ntis <- merge(ntiss.P1.anno.gencode,ntiss.P3.anno.gencode, do.normalize = TRUE)
  ntis%<>%NormalizeData()
  ntis <- merge(ntiss.P2.anno.gencode,ntis, do.normalize = TRUE)
  ntis%<>%NormalizeData()
}else if(data.use=='10Xdroplet'){
  ntis <- merge(ntiss10x.P1.anno,ntiss10x.P3.anno, do.normalize = TRUE)
  ntis%<>%NormalizeData()
  ntis <- merge(ntiss10x.P2.anno,ntis, do.normalize = TRUE)
  ntis%<>%NormalizeData()
}
Idents(ntis)=ntis@meta.data$free_annotation


ntis<-subset(ntis,idents=c(immune_clus.use$orig.clus))

genes.use<-rownames(ntis)[!grepl(pattern = "^MT-|^ERCC-|^RPS|^RPL", x = rownames(ntis))]
ntis<-ntis[genes.use,]

filter.ind<-which(immune_clus.use$orig.clus%ni%Idents(ntis))
if(length(filter.ind)>0){
  ntis<-subset(ntis,idents = immune_clus.use$orig.clus[-filter.ind])
}else{
  ntis<-subset(ntis,idents = immune_clus.use$orig.clus)
}

Idents(ntis)=factor(Idents(ntis),levels = immune_clus.use$orig.clus[immune_clus.use$orig.clus%in%unique(Idents(ntis))])
{
  sctv2 <- sctransform::vst(ntis[["RNA"]]@counts, vst.flavor = "v2")$y
  genesToUse = rownames(sctv2)
  cellsToUse = colnames(sctv2)
  ntis<-ntis[genesToUse,cellsToUse]
  all.genes <- rownames(ntis)
}

asw=c()
df<-crossing(c(200,500,1000,2000,4000),c('1hvg','2hdg','3heg'),c('singleR','tissue'))
df%<>%as.data.frame()
colnames(df)=c('nfeature','method','cluster')
df$method%<>%gsub(pattern='[0-9]',replacement='')

for(nf in c(200,500,1000,2000,4000)){
  plot.list=list()
  sil_df<-data.frame(ngene=numeric(0),
                     method=character(0),
                     clustered.by=character(0),
                     cluster=character(0),
                     cluster.size=numeric(0),
                     cluster.sil.score=numeric(0)
  )
  for(meth in c('hvg','hdg','heg')){
    if(meth=='hvg'){
      ntis<-FindVariableFeatures(ntis, selection.method = "vst", nfeatures = nf)
      
    }else if(meth=='hdg'){
      deviances <- scry::devianceFeatureSelection(as.matrix(ntis@assays$RNA@counts), fam=c("poisson"))
      o <- order(deviances, decreasing=TRUE)
      HDG <- rownames(ntis@assays$RNA@counts)[o][1:nf]
      VariableFeatures(object = ntis)=HDG
    }else if(meth=='heg'){
      cpm<-edgeR::cpm(as.matrix(ntis@assays$RNA@counts))
      genemeans<-rowMeans(cpm)
      o<-order(genemeans,decreasing = T)
      HEG<-rownames(cpm)[o][1:nf]
      Seurat::VariableFeatures(object = ntis)=HEG
    }
    
    ntis <- ScaleData(ntis, features = all.genes)# ,do.scale = F,do.center = F
    ntis <- RunPCA(ntis, features = VariableFeatures(object = ntis))
    ntis <- RunUMAP(ntis, dims = 1:NPC)
    
    ref=ref.mrna.use
    
    if(length(filter.ind)>0){
      ref<-ref[rownames(ref)%in%rownames(ntis),ref$label.refine%in%immune_clus.use$bulkseq_clus[-filter.ind]]
    }else{
      ref<-ref[rownames(ref)%in%rownames(ntis),ref$label.refine%in%immune_clus.use$bulkseq_clus]
    }
    results<-SingleR(test=as.SingleCellExperiment(ntis),ref=ref,labels=ref$label.refine,BPPARAM =  multicoreParam)
    ntis$singlr_label<-factor(results$labels,levels=immune_clus.use$bulkseq_clus[immune_clus.use$bulkseq_clus%in%unique(results$labels)])
    
    ntis@meta.data$tissue%<>%str_to_title() 
    ntis.0<-ntis
    Idents(ntis.0)=ntis$singlr_label
    if(largecluster){
      ntis.0<-subset(ntis.0,idents=table(Idents(ntis.0))[table(Idents(ntis.0))>0.01*ncol(ntis)]%>%names())#larger than 23
    }
    
    ident1<-Idents(ntis.0)
    ident2<-paste0(ntis.0@meta.data$tissue,'_',Idents(ntis.0)%>%as.vector())
    pairs_df <- expand.grid(vector1 = c(paste0(sort(unique(ntis.0@meta.data$tissue)),'_')),vector2 = levels(ident1))
    all_pairs <- apply(pairs_df, 1, function(row) paste0(row["vector1"], row["vector2"]))
    Idents(ntis.0)=factor(ident2,levels = all_pairs)
    ntis.0<-subset(ntis.0,idents=table(Idents(ntis.0))[table(Idents(ntis.0))>0.01*ncol(ntis)]%>%names())#larger than 23
    ident1.1<-Idents(ntis.0)%>%as.vector()%>%
      gsub(pattern='^Blood_',replacement='')%>%
      gsub(pattern='^Lung_',replacement='')
    Idents(ntis.0)=factor(ident1.1,levels = levels(ident1))
    
    distance_matrix <- dist(Embeddings(ntis.0[['umap']])[, 1:2])
    clusters <- as.character(Idents(ntis.0))
    c.list = names(table(clusters))
    N = length(c.list)
    for(i in 1:N){ clusters[clusters==c.list[i]] <- i }
    sil <- cluster::silhouette(clusters%>%as.numeric(), dist = distance_matrix)
    ss <-   summary(sil)
    ntis@misc$ASW <- round(ss$avg.width, digits=5)
    
    ##
    asw<-c(asw,ntis@misc$ASW)
    ##
    
    g0<-DimPlot(ntis.0, reduction = dim.reduction,cols = dark_vivid_colors[1:length(unique(Idents(ntis.0)))],
                pt.size = 0.6)
    plot.list[[paste0(meth,nf,'_by_singler')]]=g0+
      ggtitle(paste0(nf,' ', toupper(meth),'\nby singleR'),
              subtitle = paste0('silhouette width : ',round(ss$avg.width, digits=5)))+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
    
    
    sil_df.temp<-data.frame(ngene=rep(nf,N),
                            method=rep(meth,N),
                            clustered.by=rep('by singleR',N),
                            cluster=c.list,
                            cluster.size=ss$clus.sizes,
                            cluster.sil.score=ss$clus.avg.widths
    )
    sil_df<-rbind(sil_df,sil_df.temp)
    
    
    
    
    g0<-DimPlot(ntis.0,group.by = 'tissue', reduction = dim.reduction,cols = c('firebrick1','darkblue'),
                pt.size = 0.6)+theme(plot.title = element_text(hjust = 0.5))
    plot.list[[paste0(meth,nf,'_by_tissue')]]=g0+
      ggtitle(paste0(nf,' ', toupper(meth),'\nby tissue'))
    
    
    ident1<-Idents(ntis.0)
    ident2<-paste0(ntis.0@meta.data$tissue,'_',Idents(ntis.0)%>%as.vector())
    
    
    pairs_df <- expand.grid(vector1 = c(paste0(sort(unique(ntis.0@meta.data$tissue)),'_'),''),vector2 = levels(ident1))
    all_pairs <- apply(pairs_df, 1, function(row) paste0(row["vector1"], row["vector2"]))
    ap2<-apply(pairs_df%>%dplyr::filter(vector1!=''), 1, function(row) paste0(row["vector1"], row["vector2"]))
    Idents(ntis.0)=factor(ident2,levels = ap2[ap2%in%unique(ident2)])
    
    color_extend2<-color_extend[ap2%in%unique(ident2)]
    
    
    
    
    clusters <- as.character(Idents(ntis.0))
    c.list = names(table(clusters))
    N = length(c.list)
    for(i in 1:N){ clusters[clusters==c.list[i]] <- i }
    sil <- cluster::silhouette(clusters%>%as.numeric(), dist = distance_matrix)
    ss <-   summary(sil)
    ntis@misc$ASW <- round(ss$avg.width, digits=5)
    
    ##
    asw<-c(asw,ntis@misc$ASW)
    ##
    
    
    g0<-DimPlot(ntis.0, reduction = dim.reduction,cols = color_extend2[1:length(unique(Idents(ntis.0)))],
                pt.size = 0.6)+theme(plot.title = element_text(hjust = 0.5))
    plot.list[[paste0(meth,nf,'_by_tissue_singleR')]]=g0+
      ggtitle(paste0(nf,' ', toupper(meth),'\nby singleR from each tissue'),
              subtitle = paste0('silhouette width : ',round(ss$avg.width, digits=5)))+
      theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5))
    
    sil_df.temp<-data.frame(ngene=rep(nf,N),
                            method=rep(meth,N),
                            clustered.by=rep('by singleR from each tissue',N),
                            cluster=c.list,
                            cluster.size=ss$clus.sizes,
                            cluster.sil.score=ss$clus.avg.widths
    )
    sil_df<-rbind(sil_df,sil_df.temp)
  }
  
  
  title.txt = paste0(nf," features_with clusters")
  
  save(plot.list,file=paste0(dir.base,title.txt,'plot.list.RData'))
  write.table(sil_df,file=paste0(dir.base,title.txt,'_silwidth.txt'))
  ggsave(paste0(dir.base,title.txt,'.plot(by_singleR_tissue).pdf'),
         plot= grid.arrange(
           grobs = plot.list,
           widths = c(1.2, 1, 1.24),
           heights= c(1,1,1),
           layout_matrix = rbind(c(1:3),
                                 c(4:6), #NA
                                 c(7:9)),
           top = textGrob(paste0(title.txt,'(>1% total cells)'), gp = gpar(fontsize = 25))
         ),
         width = 30,height = 22
  )
}
df$asw=asw
write.table(df,file=paste0(dir.base,'total_silwidth.txt'))


# creating multiple bar plots in R
library(ggplot2)
library(magrittr)
library(RColorBrewer)

{
  main_tt='Compare_performance'
  dir.base='/custom_output_dir/'
  df<-read.table(paste0(dir.base,'total_silwidth.txt'))
  
  avg_sil_score<-df$asw
  cl<-df$cluster
  cl<-gsub(pattern='tissue',replacement = 'singleR+tissue',x=cl)
  Methods=paste0(toupper(df$method),'_',cl)
  nFeatures=df$nfeature%>%as.character()
}


# creating data frame
circle <- data.frame(avg_sil_score,Methods,nFeatures)
circle$nFeatures=factor(circle$nFeatures%>%as.character(),levels=c('200','500','1000','2000','4000'))
circle$Methods%<>%factor(levels=c('HVG_singleR','HVG_singleR+tissue','HDG_singleR','HDG_singleR+tissue','HEG_singleR','HEG_singleR+tissue'))
# # creating plot using the above data
# colorcodes=c( "#E31A1C","#1F78B4","#33A02C" , "#FB9A99","#A6CEE3" ,"#B2DF8A" )
# 
# names(colorcodes)=c('HVG_singleR','HDG_singleR','HEG_singleR','HVG_singleR+tissue','HDG_singleR+tissue','HEG_singleR+tissue')
# g1<-ggplot(circle, aes(x=nFeatures, avg_sil_score, fill = Methods)) +
#   scale_fill_manual(
#     values=colorcodes)+
#   geom_bar(stat="identity", position=position_dodge(.7),width = 0.7) +ylim(c(0,0.4))+
#   labs(title=main_tt,subtitle="Avg_silwidth")+theme(panel.background = element_rect(fill='white',colour = 'black',linewidth=1),text = element_text(size=15,color='black'))
# 
# ggsave(filename=paste0(dir.base,main_tt,'_','sil_6.pdf'),plot = g1,width=8,height=6)

colorcodes=c( "#E31A1C","#1F78B4","#33A02C" ,"#E31A1C","#1F78B4","#33A02C"  )
names(colorcodes)=c('HVG_singleR','HDG_singleR','HEG_singleR','HVG_singleR+tissue','HDG_singleR+tissue','HEG_singleR+tissue')
singleR=c('HVG_singleR','HDG_singleR','HEG_singleR')
tissue=c('HVG_singleR+tissue','HDG_singleR+tissue','HEG_singleR+tissue')



g1<-ggplot(circle%>%dplyr::filter(Methods%in%singleR), aes(x=nFeatures, avg_sil_score, fill = Methods)) +
  scale_fill_manual(
    # values=c('red','royalblue','green','#FF7276','lightblue','lightgreen'),#'darkred','darkblue','darkgreen'
    values=colorcodes)+
  geom_bar(stat="identity", position=position_dodge(.5),width = 0.5) +ylim(c(0,0.4))+
  labs(title=main_tt, subtitle = "Avg_silwidth(singleR)")+theme(panel.background = element_rect(fill='white',colour = 'black',linewidth=1),text = element_text(size=15,color='black'))

ggsave(filename=paste0(dir.base,main_tt,'_','sil_3singleR.pdf'),plot = g1,width=8,height=6)

g1<-ggplot(circle%>%dplyr::filter(Methods%in%tissue), aes(x=nFeatures, avg_sil_score, fill = Methods)) +
  scale_fill_manual(
    # values=c('red','royalblue','green','#FF7276','lightblue','lightgreen'),#'darkred','darkblue','darkgreen'
    values=colorcodes)+
  geom_bar(stat="identity", position=position_dodge(.5),width = 0.5) +ylim(c(0,0.4))+
  labs(title=main_tt, subtitle = "Avg_silwidth(singleR+tissue)")+theme(panel.background = element_rect(fill='white',colour = 'black',linewidth=1),text = element_text(size=15,color='black'))

ggsave(filename=paste0(dir.base,main_tt,'_','sil_3singleR+tissue.pdf'),plot = g1,width=8,height=6)
