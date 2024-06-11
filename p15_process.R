##Ozel p15 big result##

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

# read in------

p15_meta <- read_csv("FlyphoneDB/ozel_p15/p15_meta.csv")
celltype<-unique(p15_meta$celltype)
length(celltype)


dict_df<-read.csv("FlyphoneDB/Clusters-dict.csv")

file<-list.files(path = "FlyphoneDB/ozel_p15/output_test_dataset")
length(file)


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
        .[[p]] <= 0.01
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

#Actual process----

big_df<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(big_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")

for (f in file){
  path<-paste0("FlyphoneDB/ozel_p15/output_test_dataset/",f)
  df_o<-read.csv(path,check.names = FALSE)
  df<- data.frame(df_o[,-1], row.names = df_o[,1],check.names = FALSE)
  c_df<-make_clean_df(df)
  big_df<-rbind(big_df,c_df)
  print(f)
}

big_df<-big_df %>%
  filter(r_cluster!="character(0)")


big_df$l_cluster <- sapply(big_df$l_cluster, as.character)
big_df$r_cluster <- sapply(big_df$r_cluster, as.character)

write.csv(big_df, file = "FlyphoneDB/ozel_p15/all_p0.01.csv",row.names = FALSE)

freq<-table(p15_meta$celltype)
freq<-t(freq)
freq<-as.data.frame(freq)
freq<-freq[,-1]


#Filter lamina-----

ln_list<-c('L1','L2','L3','L4','L5','LPC')
ln_df<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(ln_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")

for(ct in ln_list){
  l_df<-big_df %>% 
    filter(grepl(ct,l_cluster))
  ln_df<-rbind(ln_df,l_df)
  
  r_df<-big_df %>% 
    filter(grepl(ct,r_cluster))
  ln_df<-rbind(ln_df,r_df)
  
}


ln_df<-ln_df %>%
  filter(l_cluster!='LLPC1')
ln_df<-ln_df %>%
  filter(r_cluster!='LLPC1')
ln_df<-ln_df %>%
  filter(l_cluster!='LPC1')
ln_df<-ln_df %>%
  filter(r_cluster!='LPC1')

ln_df[ln_df == "character(0)"]<-0
ln_df$l_cluster <- sapply(ln_df$l_cluster, as.character)
ln_df$r_cluster <- sapply(ln_df$r_cluster, as.character)



write.csv(ln_df, file = "FlyphoneDB/ozel_p15/lamina_p0.05.csv")


#Search Function-----

search_int=function(og_df,l_gene,r_gene,l_c,r_c){
  
  df<-og_df
  
  if(l_gene=='~'){
    df<-df
  }else{
    df<-df %>%
      filter(Gene_secreted == l_gene)
  }
  
  if(r_gene=='~'){
    df<-df
  }else{
    df<-df %>%
      filter(Gene_receptor == r_gene)
  }
  
  if(l_c=='~'){
    df<-df
  }else{
    df<-df %>%
      filter(l_cluster == l_c)
  }
  
  if(r_c=='~'){
    df<-df
  }else{
    df<-df %>%
      filter(r_cluster == r_c)
  }

  return(df)
}


test_df<-search_int(big_df,"~","~",'L1',"~")
