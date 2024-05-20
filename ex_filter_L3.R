# Expression filter with L3 only
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)

#main function-----------
call_l= function(x){
  return(str_split(x,'>')[[1]][1])
}
call_r= function(x){
  return(str_split(x,'>')[[1]][2])
}


express_filter<- function(pct, sd){
  r_dot_filter<-r_dot %>%
    filter(avg.exp.scaled>=sd,pct.exp>=pct)
  r_dot_filter$comb<-paste(r_dot_filter$id,r_dot_filter$features.plot,sep = '_')
  r_comb<-unique(r_dot_filter$comb)
  
  l_dot_filter<-l_dot %>%
    filter(avg.exp.scaled>=sd,pct.exp>=pct)
  l_dot_filter$comb<-paste(l_dot_filter$id,l_dot_filter$features.plot,sep = '_')
  l_comb<-unique(l_dot_filter$comb)
  
  ad5_f<-ad5 %>%
    filter (l_id %in% l_comb & r_id %in% r_comb)
  
  #print(paste0('%expression>',pct,', scaled level>',sd,', #entry=',nrow(ad5_f)))
  #return(ad5_f)
  return(nrow(ad5_f))
}

dict_df<-read.csv("FlyphoneDB/parameters/Clusters-dict.csv")
call_cluster_name = function (x){
  name<-dict_df[which(dict_df$Cluster == x),2]
  return(name)
}

LR_pair<-read.csv('FlyphoneDB/parameters/ligand_receptor_pair_new.csv')
l_list<-unique(LR_pair$Gene_secreted)
r_list<-unique(LR_pair$Gene_receptor)


############################
## iteration start point ##
###########################

#matrix to obj normalisation------

ad_obj <- readRDS("~/FlyphoneDB/ozel_p70/GSE142787_P70.rds")

ad_obj <- UpdateSeuratObject(object = ad_obj)
ad_obj<-subset(x = ad_obj, subset = L3_vs_P15 =='2') 
DefaultAssay(ad_obj) <- "RNA"
ad_obj <- NormalizeData(ad_obj)

##L3 specific------
ad_obj <- readRDS("~/FlyphoneDB/ozel_l3/GSE167266_OL.L3_P15_merged.rds")

ad_obj <- UpdateSeuratObject(object = ad_obj)
ad_obj<-subset(x = ad_obj, subset = L3_vs_P15 =='2') 
DefaultAssay(ad_obj) <- "RNA"
ad_obj <- NormalizeData(ad_obj)

l_o<-levels(ad_obj@active.ident)

##receptor--------
r_dot <- DotPlot(ad_obj,
                 features = r_list, 
                 idents = l_o)+ 
  theme(axis.text.x = element_text(angle = 90))$data
r_dot<-r_dot$data
r_dot[]<-lapply(r_dot,function(x)gsub('L3_','',x))
  
r_dot$cluster<-lapply(r_dot$id, call_cluster_name)


##ligand---------

l_dot <- DotPlot(ad_obj,
                 features = l_list, 
                 idents = l_o)+ 
  theme(axis.text.x = element_text(angle = 90))$data
l_dot<-l_dot$data
l_dot[]<-lapply(l_dot,function(x)gsub('L3_','',x))
l_dot$cluster<-lapply(l_dot$id, call_cluster_name)





#process DB result-----

ad5<-read.csv('FlyphoneDB/ozel_l3/all_p0.05.csv')

ad5$l_id<-lapply(ad5$pairing,call_l)
ad5$r_id<-lapply(ad5$pairing,call_r)


ad5$l_id<-paste(ad5$l_id,ad5$Gene_secreted,sep='_')
ad5$r_id<-paste(ad5$r_id,ad5$Gene_receptor,sep='_')



ad5_f<-express_filter(35,0.7)

stage='l3'

write.csv(ad5_f, file = (paste0("FlyphoneDB/ozel_",stage,"/wnt_fil_0.05.csv")),row.names = FALSE)
f_known<-ad5_f %>% filter(l_cluster!=''&r_cluster!='')
f_known<-f_known %>% filter(!is.na(r_cluster))
write.csv(f_known, file = (paste0("FlyphoneDB/ozel_",stage,"/wnt_fil_0.05_known.csv")),row.names = FALSE)


#sum table-----
pct_list<-c(10,20,40,60,80)
sd_list<-c(0,0.5,1,1.5,2)



df <- data.frame(matrix(ncol = 5, nrow = 5))

colnames(df) <- sd_list
rownames(df) <- pct_list

for(i in c(1:5)){
  for (j in c(1:5)){
    df[i,j]<-express_filter(pct_list[i],sd_list[j])
  }
}

