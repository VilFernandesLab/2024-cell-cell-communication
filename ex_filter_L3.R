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
r_dot[]<-lapply(r_dot,function(x)gsub('L3_','',x)) #remove 'L3' prefix
  
r_dot$cluster<-lapply(r_dot$id, call_cluster_name)


##ligand---------

l_dot <- DotPlot(ad_obj,
                 features = l_list, 
                 idents = l_o)+ 
  theme(axis.text.x = element_text(angle = 90))$data
l_dot<-l_dot$data
l_dot[]<-lapply(l_dot,function(x)gsub('L3_','',x)) #remove 'L3' prefix
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

#draw score change------

df <- data.frame(matrix(ncol = 3, nrow = 0))


for(i in c(1:5)){
  for (j in c(1:5)){
    sub_df<-express_filter(pct_list[i],sd_list[j])
    
    ad_df<-data.frame(matrix(ncol = 3, nrow = nrow(sub_df)))
    colnames(ad_df)<-c('score','percentage','sc_expression')
    ad_df$score <-sub_df$score
    ad_df$percentage <- pct_list[i] %>% as.character()
    ad_df$sc_expression <- sd_list[j] %>% as.character()
    
    df <- rbind(df, ad_df)
  }
}

ggplot(data=df, aes(x=score)) +
  #geom_histogram()+
  geom_density(adjust=1.5,color="darkblue", fill="lightblue") +
  facet_grid(vars(percentage), vars(sc_expression))+
  ggtitle('L3')+
  theme_bw()


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

