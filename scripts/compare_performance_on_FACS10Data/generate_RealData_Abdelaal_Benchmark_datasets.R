# For the Zheng sorted datasets, we downloaded the 10 PBMC-sorted populations 
# (CD14+ monocytes, CD19+ B cells, CD34+ cells, CD4+ helper T cells, CD4+/CD25+ regulatory T cells, 
#   CD4+/CD45RA+/CD25− naive T cells, CD4+/CD45RO+ memory T cells, CD56+ natural killer cells, 
#   CD8+ cytotoxic T cells, CD8+/CD45RA+ naive cytotoxic T cells) 
# from https://support.10xgenomics.com/single-cell-gene-expression/datasets; 
# https://github.com/tabdelaal/scRNAseq_Benchmark/blob/master/README.md

generate_RealData_Abdelaal_Benchmark<- function(n){
  
  ############# load real data ######################
  #dataName.list = c("Zheng_sorted")
  
  if(n==1){
    path.data = "./data/FACS10Data/"
    file.list = dir(path.data)
    
    Data = read.csv("./data/FACS10Data/Filtered_DownSampled_SortedPBMC_data.csv", row.names=1)
    Labels = as.matrix(read.csv("./data/FACS10Data/Labels.csv"))
    
    countData = as.matrix(t(Data))
    realSubtype = as.factor(Labels)
  }
  
  data = list()
  data$countData = countData
  data$realSubtype = realSubtype
  
  return(data)
}
