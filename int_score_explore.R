
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
# import data-----
ad5<-read.csv('FlyphoneDB/ozel_adult/all_p0.05.csv')
p155<-read.csv('FlyphoneDB/ozel_p15/all_p0.05.csv')

data<- readRDS('~/FlyphoneDB/ozel_adult/GSE142787_Adult.rds')

pbmc <- UpdateSeuratObject(object = data)

rm(data)

#subset to check expression-----

s_set<-subset(x = pbmc, idents = c(74,199,151,107,102,85))

AverageExpression(object = s_set,features = c('N','wry','AstA','AstA-R1','Wnt2','fz3'))$RNA

#Distribution of score-----
ggplot(ad5, aes(x=score)) + geom_histogram(binwidth=0.1)+
  scale_x_continuous(trans = 'log10')+
  scale_y_continuous(trans = 'log10')

ad5 %>%
  filter(score >= 1) %>% nrow()
# 47937

ad5 %>%
  filter(score >= 0.5) %>% nrow()
# 69397

bk<-c(1e-06,1e-05,1e-04,1e-03,1e-02,1e-01,1,10)
ggplot(ad5,aes(x=score)) + geom_step(aes(y=after_stat(y)),stat="ecdf")+
  scale_x_continuous(breaks=bk,trans = 'log10')+
  scale_y_continuous(breaks=seq(0,1,by=0.05))

ggplot(p155,aes(x=score)) + geom_step(aes(y=after_stat(y)),stat="ecdf")+
  scale_x_continuous(breaks=bk,trans = 'log10')+
  scale_y_continuous(breaks=seq(0,1,by=0.05))



#matrix to obj normalisation------
## import data-------
ad_obj<-readRDS('~/FlyphoneDB/ozel_adult/GSE142787_Adult.rds')

ad_obj <- UpdateSeuratObject(object = ad_obj)
DefaultAssay(ad_obj) <- "RNA"
ad_obj <- NormalizeData(ad_obj)

lamina<-c(108,110,109,107,117)
glia<-c(186:189,198:208)
RidgePlot(ad_obj,features = 'repo', idents = glia)

RidgePlot(ad_obj,features = c('wry','N'), idents = glia)
p<-RidgePlot(ad_obj,features = c('wry','N'), idents = c(186:189,198:208,29))
p+geom_density_ridges(bandwidth = 0.1)

DotPlot(ad_obj,features = c('wry','N'), idents = glia)
VlnPlot(ad_obj,features = c('wry','N'), idents = glia)

DotPlot(ad_obj,features = 'ser', idents = c(1:10))
VlnPlot(ad_obj,features = c('wry'), idents = c(1:20))


#violin plot -----
#Ser not found

lr_pair<-read.csv("FlyphoneDB/parameters/ligand_receptor_pair_new.csv")

id_list=list()
for (i in c(0:12)){
  id_list<-c(id_list,list((i*20):(i*20+19)))
}


##repeat----------- 

plot_v<-function(gene){
  plot_list = list()
  
  
  for (i in c(1:13)) {
    p = VlnPlot(ad_obj,features = gene, idents = id_list[[i]])
    plot_list[[i]] = p
  }
  
  pdf(paste0('~/FlyphoneDB/ozel_adult/plot/v_',gene,'.pdf'), width=10, height=5)
  print(plot_list)
  dev.off()
}





plot_v('dally')

egfr<- c("aos","grk", "Krn", "spi","vn",'boss',"Egfr",'sev')


gene<-'sev'
plot_list = list()


for (i in c(1:13)) {
  p = VlnPlot(ad_obj,features = gene, idents = id_list[[i]])
  plot_list[[i]] = p
}

pdf(paste0('~/FlyphoneDB/ozel_adult/plot/v_',gene,'.pdf'), width=10, height=5)
plot_list
dev.off()
plot_list = list()


for (i in c(1:13)) {
  p = VlnPlot(ad_obj,features = gene, idents = id_list[[i]])
  plot_list[[i]] = p
}

pdf(paste0('~/FlyphoneDB/ozel_adult/v_',gene,'.pdf'), width=10, height=5)
plot_list
dev.off()

#Dotplot-----
Pathway_core_components <- read.csv("FlyphoneDB/parameters/fly_phone_pathway_core_components.csv")
##EGFR------

e_df <- subset(Pathway_core_components, Pathway == 'EGFR signaling pathway')

egfr<- c("grk", "Krn", "spi","vn")


FeaturePlot(ad_obj,features = 'N',max.cutoff = 1,label = TRUE)


#big dot plot------
LR_pair<-read.csv('FlyphoneDB/parameters/ligand_receptor_pair_new.csv')
l_list<-unique(LR_pair$Gene_secreted)
r_list<-unique(LR_pair$Gene_receptor)
DotPlot(ad_obj,features = l_list, idents = c(1:40))+ theme(axis.text.x = element_text(angle = 90))


id_list_d=list()
for (i in c(0:5)){
  id_list_d<-c(id_list_d,list((i*40):(i*40+39)))
}

plot_list = list()


for (i in c(1:6)){
  p = DotPlot(ad_obj,features = r_list, idents = id_list_d[[i]])+ theme(axis.text.x = element_text(angle = 90))
  plot_list[[i]] = p
  
}

pdf(paste0('~/FlyphoneDB/ozel_adult/plot/d_','receptor','.pdf'), width=15, height=10)
plot_list
dev.off()


## extract data ------
call_cluster_name = function (x){
  name<-dict_df[which(dict_df$Cluster == x),2]
  return(name)
}
dict_df<-read.csv("FlyphoneDB/parameters/Clusters-dict.csv")
##receptor--------
r_dot <- DotPlot(ad_obj,
                     features = r_list, 
                     idents = c(1:235))+ 
  theme(axis.text.x = element_text(angle = 90))$data
r_dot<-r_dot$data
r_dot$cluster<-lapply(r_dot$id, call_cluster_name)


r_dot_filter<-r_dot %>%
  filter(avg.exp.scaled>=1,pct.exp>=40)
r_dot_filter$comb<-paste(r_dot_filter$id,r_dot_filter$features.plot,sep = '_')
r_comb<-unique(r_dot_filter$comb)
##ligand---------

l_dot <- DotPlot(ad_obj,
                 features = l_list, 
                 idents = c(1:235))+ 
  theme(axis.text.x = element_text(angle = 90))$data
l_dot<-l_dot$data
l_dot$cluster<-lapply(l_dot$id, call_cluster_name)


l_dot_filter<-l_dot %>%
  filter(avg.exp.scaled>=1,pct.exp>=40)
l_dot_filter$comb<-paste(l_dot_filter$id,l_dot_filter$features.plot,sep = '_')
l_comb<-unique(l_dot_filter$comb)


##expression filter-----

call_l= function(x){
  return(str_split(x,'>')[[1]][1])
}
call_r= function(x){
  return(str_split(x,'>')[[1]][2])
}


ad5<-read.csv('FlyphoneDB/ozel_adult/p0.05_spatial_fil.csv')
ad5<-read.csv('FlyphoneDB/ozel_adult/all_p0.05.csv')

ad5$l_id<-lapply(ad5$pairing,call_l)
ad5$r_id<-lapply(ad5$pairing,call_r)


ad5$l_id<-paste(ad5$l_id,ad5$Gene_secreted,sep='_')
ad5$r_id<-paste(ad5$r_id,ad5$Gene_receptor,sep='_')



ad5_f<-ad5 %>%
  filter (l_id %in% l_comb & r_id %in% r_comb)


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
  return(ad5_f)
}

ad5_f<-express_filter(80,2)
write.csv(ad5_f, file = ("FlyphoneDB/ozel_adult/all_80_2.csv"),row.names = FALSE)

pct_list<-c(10,20,40,60,80)
sd_list<-c(0,0.5,1,1.5,2)

for(i in c(1:5)){
  for (j in c(1:5)){
    count_filter(pct_list[i],sd_list[j])
  }
}

