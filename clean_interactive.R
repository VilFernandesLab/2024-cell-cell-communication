##Cluster correspondence ##


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

# input df -------

## dictionary df -------
dict_df<-read.csv("FlyphoneDB/Clusters-dict.csv")
dict_df<-dict_df[c(1,2,183),]
## intreactive scores df -------
df<-read.csv("FlyphoneDB/1,3,198_subset/interaction_list.csv",check.names = FALSE)
LR_pairs_one<- data.frame(df[,-1], row.names = df[,1],check.names = FALSE)


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
        .[[p]] <= 0.05
      ) 
    
    fil_df<-fil_df[,c(1:3,s)]
    colnames(fil_df)[ncol(fil_df)]<-'score'
    fil_df['pairing']<-str_split(colnames(df)[s],'_')[[1]][1]
    
    clean_df<-rbind(clean_df,fil_df)
  }
  
  
  
  clean_df$l_cluster<-lapply(clean_df$pairing,call_l)
  clean_df$l_cluster<-lapply(clean_df$l_cluster, call_cluster_name)
  
  clean_df$r_cluster<-lapply(clean_df$pairing,call_r)
  clean_df$r_cluster<-lapply(clean_df$r_cluster, call_cluster_name)
  
  return(clean_df)
}

test_clean<-make_clean_df(LR_pairs_one)
