#4 HVGs ARI in 1 page - filter rate ARI plot

#ARI plot

#source("~/3.scripts/main_figures/supple_fig.9_filter_rate.R")
#generate_supple_fig.9_filter_rate(runName = "DrHiDepth_featureNumbers", dataName.list= "alz_scd_v3")

generate_supple_fig.9_filter_rate <- function( runName = "DrHiDepth_featureNumbers", dataName.list= "alz_scd_v3"){
  
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
  dataName.list =  "alz_scd_v3"
  
  #1. data 
  N.data = length(dataName.list)
  #2. dataName 
  for(n in 1:N.data){
    n = 1
    print(paste("n:",n))
    
    #3. filter_rate
    filter.list = c(0, 0.2, 0.5, 0.7) #c(0, 0.1, 0.2) #, 0.5, 0.7)
    N.f = length(filter.list)
    rm(filter.table)
    for(f in 1:N.f){
      print(paste("f:",f))
      filter_rate = filter.list[f]
      
      tableName = paste("DrHiDepth",dataName.list[n], filter_rate, sep="_")
      ARIfileName = paste(file.path.result, "ARI.table", tableName,".txt", sep="")
      
      models.list = c("scTv2") #c( "scTv2","logN")
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
          #method.list = c("scTv2.HVG.vst", "scTv2.HVG.disp", "scTv2.HVG.plus.HDG", "scTv2.HVG.plus.HEG" )
          #method.list = c("scTv2.HVG.vst", "scTv2.HVG.plus.HDG", "scTv2.HVG.plus.HEG" )
          method.list = c("scTv2.HVG.vst", "scTv2.HDG", "scTv2.HEG.cpm" )
          
          
          w = which(m %in% method.list )
          #w <- grep(paste(models.list[k],sep=""), m)
          #print(paste("k:", k, " - model:", models.list[k]))
          score.table = tmp[,w]
          colnames(score.table)=paste(colnames(score.table),".filter" ,filter_rate, sep="")
          
          if(f==1){ filter.table = score.table
          }else if(f>1){ filter.table = cbind(filter.table, score.table)}
          
          # # score plots
          # plot.table = as.matrix(score.table)
          # method.list  = rep(rownames(plot.table), ncol(plot.table))
          # score.list = c(plot.table) #sprintf("%.3f",c(plot.table))
          # 
          # plot.data = data.frame(method=method.list, nfeatures=nfeatures.list, ARI=score.list)
          # plot.data$nfeatures = factor(plot.data$nfeatures, levels=rownames(plot.table))
          # ptitle = models.list[k]
          # 
          # p <- ggplot(plot.data,  aes(x=nfeatures, y=ARI, color=method, group=method)) + 
          #   geom_point(col=3) + 
          #   geom_line() + 
          #   ggtitle(ptitle) + theme_bw() + 
          #   theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=15),
          #         axis.title = element_text(size=15, face="bold"),
          #         plot.title = element_text(size=15, hjust=0.5, face="bold"),
          #         legend.title = element_text(size=14),
          #         legend.text=element_text(size=12))
          # 
          # pName = paste("p",s,k, sep="")
          # print(paste("p.s.k:",pName,sep=""))
          # assign(pName, p)
        }#s1
      }#k for model 
      
      
    }#f for filter
    
    filter.table = filter.table[,c(1,2,3,4,7,10)]
    
    # score plots
    score.table = filter.table
    plot.table = as.matrix(score.table)
    method.list = rep(gsub(paste(models.list[k],".",sep=""), "",colnames(plot.table)), each=nrow(plot.table))
    nfeatures.list = rep(rownames(plot.table), ncol(plot.table))
    score.list = c(plot.table) #sprintf("%.3f",c(plot.table))
    
    plot.data = data.frame(method=method.list, nfeatures=nfeatures.list, ARI=score.list)
    plot.data$nfeatures = factor(plot.data$nfeatures, levels=rownames(plot.table))
    #l=c("HDG", "HEG.cpm", "HVG.vst")
    #f=paste(models.list[k], l, sep=".")
    #plot.data$method = factor(plot.data$method, levels=l)
    
    ptitle = models.list[k]

    p <- ggplot(plot.data,  aes(x=nfeatures, y=ARI, color=method, group=method, linetype=method, size=method)) +
      #geom_point(col=1) +
      geom_line() +
      ggtitle(ptitle) + theme_bw() +
      theme(axis.text=element_text(angle=45, hjust=1, vjust=1, color="black", size=14),
            axis.title = element_text(size=15, face="bold"),
            plot.title = element_text(size=15, hjust=0.5, face="bold"),
            legend.title = element_text(size=14),
            legend.text=element_text(size=12))
    
    
    m = unique(plot.data$method)
    #display.brewer.all()
    # user_palette[grep("HDG",m)]<- brewer.pal(9,"Set1")[1]
    # user_palette[grep("HEG.cpm",m)]<- brewer.pal(9,"Set1")[8] #4 pupple
    # user_palette[grep("HVG.vst",m)]<- brewer.pal(9,"Blues")[8]
    
    user_palette = c()
    user_palette[grep("HVG.vst.filter0",m)]<- brewer.pal(9,"Blues")[8]
    user_palette[grep("HVG.vst.filter0.2",m)]<- brewer.pal(9,"Blues")[7]
    user_palette[grep("HVG.vst.filter0.5",m)]<- brewer.pal(9,"Blues")[6]
    user_palette[grep("HVG.vst.filter0.7",m)]<- brewer.pal(9,"Blues")[5]
    
    user_palette[grep("HDG.filter0",m)]<-brewer.pal(9,"Set1")[1]
    user_palette[grep("HEG.cpm.filter0",m)]<- brewer.pal(9,"Set1")[8] 
   
    plot_col= user_palette
    plot_lty= rep('solid',length(m))
    plot_size=rep(2,length(m))

    names(plot_col) = m
    names(plot_lty)=names(plot_size)=names(plot_col)

    plot_lty[['HVG.vst.filter0.7']]="4121"

    plot_size[['HVG.vst.filter0.2']]=1.5
    plot_size[['HVG.vst.filter0.5']]=0.9
    plot_size[['HVG.vst.filter0.7']]=0.7

    #st=plot_style(m=unique(plot.data$method))
    
    p = p + scale_linetype_manual(values = plot_lty)+
      scale_colour_manual(values = plot_col)+
      scale_size_manual(values = plot_size)+
      guides(guide_legend(nrow = 3, byrow = TRUE))
    

    
    #******* summary.result plot
    # integrated ARI DEG plot
    file.path.summary = "/store/juok/research/Research__scNormalization/results/summary.plot/"
    timestamp = paste0("___",format(Sys.time(), "%Y%m%d_%H%M%S"))
    plotName = paste("Supple_Fig9.","ARI_in_diff_filter_rate_", dataName.list[n], "_",timestamp,".pdf", sep="")
    
    pdfName = paste(file.path.summary, plotName, sep="")
    
    default.txt = paste("ARI in different filter rates in alz scd data")
    
    pdf(pdfName, width=10, height=10)
    
    library(gridExtra)
    library(grid)
    title.txt = default.txt 
    # pskd
    grid.arrange(p,
                 top = textGrob(title.txt, gp = gpar(fontsize = 25)),
                 nrow = 1)
    
    dev.off()
    message(paste("generated summary plot!..\n",pdfName))
  }#n for dataN
  
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

