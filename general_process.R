# general result process#

# library -------
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(future.apply)
  library(Seurat)
  library(RColorBrewer)
  library(reshape2)
  library(network)
  library(igraph)
  library(dplyr)
})




#functions----------
call_l= function(x){
  return(str_split(x,'>')[[1]][1])
}
call_r= function(x){
  return(str_split(x,'>')[[1]][2])
}

call_cluster_name = function (x){
  name<-dict_df[which(dict_df$Cluster == x),2]
  return(name)
}


make_clean_df = function (df){
  clu_sqr=(ncol(df)-3)/2
  clean_df<-data.frame(matrix(ncol = 5, nrow = 0))
  colnames(clean_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")
  for (x in c(1:clu_sqr)) {
    s=(x-1)*2+4
    p=(x-1)*2+5
    
    fil_df<- df%>%
      filter(
        .[[p]] <= 0.1
      ) 
    
    if (nrow(fil_df)==0){
      
    }else{
      fil_df<-fil_df[,c(1:3,s)]
      colnames(fil_df)[ncol(fil_df)]<-'score'
      fil_df['pairing']<-str_split(colnames(df)[s],'_')[[1]][1]
      
      clean_df<-rbind(clean_df,fil_df)
    }
  }
  
  clean_df$l_cluster<-lapply(clean_df$pairing,call_l)
  clean_df$l_cluster<-lapply(clean_df$l_cluster, call_cluster_name)
  
  clean_df$r_cluster<-lapply(clean_df$pairing,call_r)
  clean_df$r_cluster<-lapply(clean_df$r_cluster, call_cluster_name)
  
  return(clean_df)
}

# read in------
stage='l3'
meta <- read.csv(paste0("FlyphoneDB/ozel_",stage,"/meta.csv"))
celltype<-unique(meta$celltype)
length(celltype)


dict_df<-read.csv("FlyphoneDB/parameters/Clusters-dict.csv")
dict_df$Assignment <- gsub("*","",as.character(dict_df$Assignment),fixed = TRUE)
file<-list.files(path = paste0("FlyphoneDB/ozel_",stage,"/output_test_dataset"))
length(file)

#Actual process----

big_df<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(big_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")

for (f in file){
  path<-paste0("FlyphoneDB/ozel_",stage,"/output_test_dataset/",f)
  df_o<-read.csv(path,check.names = FALSE)
  df<- data.frame(df_o[,-1], row.names = df_o[,1],check.names = FALSE)
  c_df<-make_clean_df(df)
  big_df<-rbind(big_df,c_df)
  print(f)
}

big_df<-big_df %>%
  filter(r_cluster!="character(0)")
big_df<-big_df %>%
  filter(l_cluster!="character(0)")


big_df$l_cluster <- sapply(big_df$l_cluster, as.character)
big_df$r_cluster <- sapply(big_df$r_cluster, as.character)

write.csv(big_df, file = paste0("FlyphoneDB/ozel_",stage,"/all_p0.05.csv"),row.names = FALSE)



#read in (for entry length)------
stage<-'p15'
df<-read.csv(paste0("FlyphoneDB/ozel_",stage,"/all_p0.05.csv"))
nrow(df)
