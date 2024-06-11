#Search big table----
library(dplyr)
library(tidyverse)

#discarded cells have the identity ‘0’


##Import data-----

p15_1<-read.csv("FlyphoneDB/ozel_p15/lamina_p0.01.csv")

p15_5<-read.csv("FlyphoneDB/ozel_p15/lamina_p0.05.csv")#Insert your file path here
p15_5<-p15_5[,2:8]
ad_1<-read.csv("FlyphoneDB/ozel_adult/lamina_p0.01.csv")
ad_1<-ad_1[,3:9]
ad_5<-read.csv("FlyphoneDB/ozel_adult/lamina_p0.05.csv")
ad_5<-ad_5[,3:9]
###read in all----
ad<-read.csv("FlyphoneDB/ozel_adult/all_p0.05.csv")
p70<-read.csv("FlyphoneDB/ozel_p70/all_p0.05.csv")
p50<-read.csv("FlyphoneDB/ozel_p50/all_p0.05.csv")
p40<-read.csv("FlyphoneDB/ozel_p40/all_p0.05.csv")
p30<-read.csv("FlyphoneDB/ozel_p30/all_p0.05.csv")
p15<-read.csv("FlyphoneDB/ozel_p15/all_p0.05.csv")
l3<-read.csv("FlyphoneDB/ozel_l3/all_p0.05.csv")

dict_df<-read.csv("FlyphoneDB/parameters/Clusters-dict.csv")
dict_df$Assignment <- gsub("*","",as.character(dict_df$Assignment),fixed = TRUE)

p15_gl<-read.csv("FlyphoneDB/ozel_p15/p15_glia_lamina.csv")
##Function-----

call_cluster_name = function (x){
  name<-dict_df[which(dict_df$Cluster == x),2]
  return(name)
}

call_l= function(x){
  return(str_split(x,'>')[[1]][1])
}
call_r= function(x){
  return(str_split(x,'>')[[1]][2])
}

search_int=function(og_df,l_gene,r_gene,l_c,r_c){ 
  
  #og_df: the dataframe you want to filter
  #l_gene: Ligend gene you are looking for
  #r_gene: Receptor gene you are looking for
  #l_c : cell type where the ligend gene is expressed
  #r_c : cell type where the receptor gene is expressed
  # IF you don't have a specific parameter you are looking for,
  # PLEASE put '~'
  
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

#EXAMPLE-----

test_df<-unique(search_int(p15_5,"~","N",'L1',"~") )
test_2<-unique(search_int(ad_5,"~","N",'L1',"~") )
#save the filtered dataframe to a new data frame so you can view it

# Glia----
f_list<-c('Astrocyte','Epithelial','cartridge','Ensheathing','Cortex','Satellite','Chiasm',
          'Subperi','Chalice','Perineurial')

#Adult (P<0.05)------

ad_5$l_cluster<-lapply(ad_5$pairing,call_l)
ad_5$l_cluster<-lapply(ad_5$l_cluster, call_cluster_name)

ad_5$r_cluster<-lapply(ad_5$pairing,call_r)
ad_5$r_cluster<-lapply(ad_5$r_cluster, call_cluster_name)

ad_5$l_cluster <- sapply(ad_5$l_cluster, as.character)
ad_5$r_cluster <- sapply(ad_5$r_cluster, as.character)

write.csv(ad_5, file = "FlyphoneDB/ozel_adult/lamina_p0.05.csv")

g_df<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(g_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")

for(ct in f_list){
  l_df<-ad_5 %>% 
    filter(grepl(ct,l_cluster))
  g_df<-rbind(g_df,l_df)
  
  r_df<-ad_5 %>% 
    filter(grepl(ct,r_cluster))
  g_df<-rbind(g_df,r_df)
  
}

g_df$l_cluster <- sapply(g_df$l_cluster, as.character)
g_df$r_cluster <- sapply(g_df$r_cluster, as.character)
write.csv(g_df, file = "FlyphoneDB/ozel_adult/ad5_glia_lamina.csv")

g_df<-read.csv("FlyphoneDB/ozel_adult/ad5_glia_lamina.csv")

#P15 (P<0.05)------

p15_5$l_cluster<-lapply(p15_5$pairing,call_l)
p15_5$l_cluster<-lapply(p15_5$l_cluster, call_cluster_name)

p15_5$r_cluster<-lapply(p15_5$pairing,call_r)
p15_5$r_cluster<-lapply(p15_5$r_cluster, call_cluster_name)

g15_df<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(g15_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")

for(ct in f_list){
  l_df<-p15_5 %>% 
    filter(grepl(ct,l_cluster))
  g15_df<-rbind(g15_df,l_df)
  
  r_df<-p15_5 %>% 
    filter(grepl(ct,r_cluster))
  g15_df<-rbind(g15_df,r_df)
  
}

g15_df$l_cluster <- sapply(g15_df$l_cluster, as.character)
g15_df$r_cluster <- sapply(g15_df$r_cluster, as.character)
write.csv(g15_df, file = "FlyphoneDB/ozel_p15/p15_glia_lamina.csv")


# Astrocyte------

as_l<-ad_5 %>% 
  filter(grepl('Astrocyte',l_cluster))


as_70<-p70 %>% 
  filter(grepl('Astrocyte',l_cluster) | grepl('Astrocyte',r_cluster))
write.csv(as_70, file = "FlyphoneDB/ozel_p70/p70_int_Astrocyte.csv")


as_50<-p50 %>% 
  filter(grepl('Astrocyte',l_cluster) | grepl('Astrocyte',r_cluster))
write.csv(as_50, file = "FlyphoneDB/ozel_p50/p50_int_Astrocyte.csv")

#adult 0.1------
ad_1$l_cluster<-lapply(ad_1$pairing,call_l)
ad_1$l_cluster<-lapply(ad_1$l_cluster, call_cluster_name)

ad_1$r_cluster<-lapply(ad_1$pairing,call_r)
ad_1$r_cluster<-lapply(ad_1$r_cluster, call_cluster_name)

ad_1$l_cluster <- sapply(ad_1$l_cluster, as.character)
ad_1$r_cluster <- sapply(ad_1$r_cluster, as.character)

write.csv(ad_1, file = "FlyphoneDB/ozel_adult/lamina_p0.01.csv")