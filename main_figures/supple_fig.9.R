# sfig.9.a.b.c.d (DiHiData- alz scd) ARI,DEG * scTv2, logN- filter 10% 20% 50% 70%

#source("~/3.scripts/main_figures/supple_fig.9.R")
#generate_supple_fig.9(runName = "DrHiDepth_featureNumbers.hvg.norm", dataName.list= "alz_scd_v3")

generate_supple_fig.9 <- function( runName = "DrHiDepth_featureNumbers", dataName.list= "alz_scd_v3"){
  library("ggplot2")
  library(RColorBrewer)
  library(scales)
  
  dir.result = "tmp.sim.DrHiDepth_featureNumbers"
  runName="DrHiDepth_featureNumbers"
  
  file.path.result = paste("/store/juok/research/Research__scNormalization/results/", dir.result,  "/",sep = "")
  
  
  #DrHiDepthSimData v3&v4
  path.data = "/store/juok/research/Research__scNormalization/data/Dr.Hi.simData/"
  file.list = dir(path.data) 
  w = grep(".data.", file.list)
  datafile.list = file.list[w]
  dataNameAll.list = gsub(".data.RDS", "", datafile.list)
  w3=grep("_v3", dataNameAll.list)
  dataName_v3.list = dataNameAll.list[w3]
  w4=grep("_v4", dataNameAll.list)
  dataName_v4.list = dataNameAll.list[w4]
  
  #1. data v3 /v4
  version.data.list = list(dataName_v3.list) #, dataName_v4.list)
  N.v = length(version.data.list)
  for(v in 1:N.v){
    print(paste("v:", v))
    #dataName.list =version.data.list[[v]]
    dataName.list =  "alz_scd_v3"
    N.data = length(dataName.list)
    
    #2. dataName 
    for(n in 1:N.data){
      #n = 1
      print(paste("n:",n))
      
      #3. filter_rate
      filter.list = c (0.1, 0.2, 0.5, 0.7) #c(0) 
      N.f = length(filter.list)
      for(f in 1:N.f){
        print(paste("f:",f))
        filter_rate = filter.list[f]
        
        tableName = paste("DrHiDepth",dataName.list[n], filter_rate, sep="_")
        ARIfileName = paste(file.path.result, "ARI.table", tableName,".txt", sep="")
        DEGallfileName = paste(file.path.result, "DEGall.table", tableName,".txt", sep="")
        DEGupfileName = paste(file.path.result, "DEGup.table", tableName,".txt", sep="")
        DEGdownfileName = paste(file.path.result, "DEGdown.table", tableName,".txt", sep="")
        EGRG1fileName = paste(file.path.result, "EGR.G1.table", tableName,".txt", sep="")
        EGRG2fileName = paste(file.path.result, "EGR.G2.table", tableName,".txt", sep="")
        EGRG3fileName = paste(file.path.result, "EGR.G3.table", tableName,".txt", sep="")
        EGRG4fileName = paste(file.path.result, "EGR.G4.table", tableName,".txt", sep="")
        EGRG5fileName = paste(file.path.result, "EGR.G5.table", tableName,".txt", sep="")
        
        models.list = c("scTv2","logN")
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
                geom_line() + coord_cartesian(ylim = c(0, 0.8)) +
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
              dfileName.list = c(DEGallfileName, DEGupfileName, DEGdownfileName)
              dytitle.list = c("All DEG True Rate", "Up DEG True Rate", "Down DEG True Rate")
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
                  geom_line() + coord_cartesian(ylim = c(0.47, 0.65)) +
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
        
        #******* summary.result plot
        # integrated ARI DEG plot
        filter.txt = 100*(1-filter_rate)
        file.path.summary = "/store/juok/research/Research__scNormalization/results/summary.plot/"
        timestamp = paste0("___",format(Sys.time(), "%Y%m%d_%H%M%S"))
        plotName = paste("Supple_Fig.9",letters[f],"__FeatureSelectionTrend_in_ARI_&_DEG","_DrHiDepthSimData_", dataName.list[n], "_", "Hexp",filter.txt,"pro" ,"_",timestamp,".pdf", sep="")
        
        pdfName = paste(file.path.summary, plotName, sep="")
        
        default.txt = paste("Feature Selection Trend in ARI & DEG",  "\n",
                            "DrHiSimData_", dataName.list[n], "HighlyExp", filter.txt,"%", sep="")
        
        pdf(pdfName, width=20, height=15)
        
        library(gridExtra)
        library(grid)
        title.txt = default.txt 
        # ps(style ARI/DGR)k(model-scTv2/logN)d(deg all/up/down)
        grid.arrange(p11, p12,  
                     p211, p221, 
                     top = textGrob(title.txt, gp = gpar(fontsize = 25)),
                     nrow = 2)
        
        dev.off()
        message(paste("generated summary plot!..\n",pdfName))
      }#f for filter
    }#n for dataN
  }#v for version data
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

