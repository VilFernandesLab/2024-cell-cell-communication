# NON essential
# additional analysis for the number of entries in each signalling pathway.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)


#ligand_receptor pair frequency-----
score$gene_pair<-paste(score$Gene_secreted,score$Gene_receptor,sep = '-')
pair_frq <- score %>%
  group_by(pathway_receptor, gene_pair) %>%
  summarize(count=n())

# visualisation
ggplot(pair_frq, aes(x="", y=count, fill=pathway_receptor)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)


pair_frq <- score %>%
  group_by(pathway_receptor, gene_pair) %>%
  summarize(count=n())


#pathway frequency----

p5<-read.csv('FlyphoneDB/ozel_p15/all_p0.05.csv') 

path_frq <- p5 %>%
  group_by(pathway_receptor, pairing) %>%
  summarize(count=n())
table(ad5$pathway_receptor)
## functions------
dict_df<-read.csv("FlyphoneDB/parameters/Clusters-dict.csv")

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


##rearrange------
# give you a nicer and more inform data frame! Try it on!
path_frq$l_cluster<-lapply(path_frq$pairing,call_l)
path_frq$l_cluster<-lapply(path_frq$l_cluster, call_cluster_name)

path_frq$r_cluster<-lapply(path_frq$pairing,call_r)
path_frq$r_cluster<-lapply(path_frq$r_cluster, call_cluster_name)


path_frq<-path_frq[,c(1,4,5,3)]

path_frq$l_cluster <- sapply(path_frq$l_cluster, as.character)
path_frq$r_cluster <- sapply(path_frq$r_cluster, as.character)
