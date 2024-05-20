##hierarchical clustering##

#Load Library------
library(Seurat)
library(SeuratObject)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(dendextend)
library(reshape2)
library(ggdendro)
library(grid)


#Adult Glia-LN------

ad_obj <- readRDS("~/FlyphoneDB/ozel_adult/GSE142787_Adult.rds") 
ad_obj <- UpdateSeuratObject(object = ad_obj)
##L3 separation-----
pbmc <- readRDS("~/FlyphoneDB/ozel_l3/GSE167266_OL.L3_P15_merged.rds")#read in RDS
pbmc <- UpdateSeuratObject(object = pbmc)
ad_obj<-subset(x = pbmc, subset = L3_vs_P15 =='2') 


DefaultAssay(ad_obj) <- "RNA"


aver.exp <- AverageExpression(ad_obj) #Average expression
aver.exp.frame <- as.data.frame(aver.exp[['RNA']])



mydata <-aver.exp.frame
mydata_transposed <- as.data.frame(t(as.matrix(mydata)))
mydata_scale <- as.data.frame(scale(mydata_transposed))
mydata_dist <- dist(mydata_scale, method = 'euclidean') #hier clustering with euclidean score
# calculate hierarchical clustering and plot
clust_mydata <- hclust(mydata_dist, method = "complete") %>% as.dendrogram()


clust_mydata <- set(clust_mydata , "labels_cex", 0.5)
plot(clust_mydata)

# check neuron/gila------


g_list<-labels(clust_mydata)
g_list<-gsub('g','',g_list) %>% as.numeric() #order of cluster as a list
#not necessary!!!!


#pull list of average expression of gene markers
f_df<-aver.exp.frame[c('elav','repo','Rdl','Frq1','Nckx30C','fne','CG32032', 'AnxB9','GstE12','onecut','dpn','shg'),]
#save order/level of the cluster from dendrogram
o_l<-order.dendrogram(clust_mydata)
#reorder f_df in hierarchical order
f_o<-f_df[,o_l]

#write.csv(f_o, file = ("FlyphoneDB/ozel_adult/hie_gene_pre.csv"),row.names = TRUE)



##heatmap--------

#modify the average expression data frame for heat map plotting
f_m<-data.matrix(f_o)
f_long<-melt(f_m)
colnames(f_long) <- c("x", "y", "value")
#Plot the heat map
p1<-ggplot(f_long, aes(x = y, y = x, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(20, "Blue-Red"),
                       trans = 'log' ) +
  #geom_vline(xintercept = 'g3',hjust=1)+ 
  #straight line for putative glia/neuron boundary, 
  #please adjust based on your own result
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  labs(x='Cluster in Hierachical order', y="gene")


ddata_x <- dendro_data(clust_mydata )

theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank()
  #axis.ticks.length = element_blank()
)

p2 <- ggplot(segment(ddata_x)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  theme_none + theme(axis.title.x=element_blank())

 
grid.newpage()
print(p1, vp=viewport(0.8, 0.8, x=0.4, y=0.4))
print(p2, vp=viewport(0.8, 0.2, x=0.4, y=0.9))

##sum_df-------

# make summary DF with hierarchical clustering order
f_o<-t(f_o)
f_o<- cbind(cluster = rownames(f_o), f_o)
f_o<-as.data.frame(f_o)
f_o$cluster<-gsub('g','',f_o$cluster)%>% as.numeric()
f_o<-f_o[c(2:200),] #this is to remove the pseudo cartilage might be adult data set specific

#cluster to cell type correspondance
dict_df<-read.csv("FlyphoneDB/parameters/Clusters-dict.csv")
call_cluster_name = function (x){
  name<-dict_df[which(dict_df$Cluster == x),2]
  return(name)
}
f_o$id<-lapply(f_o$cluster,call_cluster_name)
#reset index number
row.names(f_o)<-c(1:199)

#assign putative neuron and glia base on hierarchical clustering
f_o$hier<-'neuron' #first assume all nueron
f_o[c(1:18),'hier']<-'glia' #later assign glia, change index number for your own result
f_o<-f_o[,c(1,13,14,2:12)] #re-order column
f_o$id<-f_o$id %>%as.character()
write.csv(f_o, file = ("FlyphoneDB/ozel_adult/hie_gene.csv"),row.names = TRUE) #write out 

#WNT4/10------
#for Ramya, please ignore this
w_dot<-DotPlot(ad_obj,
                 features = c('Wnt4','Wnt10'), 
                 idents = c(227,228))
w_data<-w_dot$data

w_obj<-subset(ad_obj,idents = c(227,228))
w_obj<-ScaleData(w_obj, features =c('Wnt4','Wnt10') )
DoHeatmap(object = w_obj,features = c('Wnt4','Wnt10'))
w_dot
