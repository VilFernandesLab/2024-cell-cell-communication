# Original script------

# open packages
library(tidyr)
library(dplyr)
library(ggtext)
library(magrittr)
# open file
data <- readRDS('~/FlyphoneDB/1,3,198_subset/1,3,198.rds')

Pathway_core_components <- read.csv("FlyphoneDB/parameters/fly_phone_pathway_core_components.csv")

unique(Pathway_core_components$Pathway)
df <- subset(Pathway_core_components, Pathway == 'NOTCH signaling pathway')
genes <- df$Gene.Symbol
# indication to use gene expression on the RNA assay
DefaultAssay(data) <- "RNA"
# plot DotPlot with 2 y axis (left axis with cluster numbers and right axis with corresponding glial subtype names)
# run the DotPlot for marker genes
plot <- DotPlot(data, scale= T,features = c('CG6126','CG9743','Tret1-1',
                                            'JhI-21','Cyp311a1','CG7135','Fas3','CG13003',
                                            'GstE9','ltl','NimB4','DAT','wrapper','Rfabg',
                                            'Dop1R1','ana','Ork1','Tsf1','Obp18a','Dtg','GILT1','Ndae1',
                                            'CG43795','axo','List','GstT4',
                                            'htl','mbc','Tet','wun2','Act79B'))


plot <- DotPlot(data, scale= T,features = genes)+ theme(axis.text.x = element_text(angle = 90))
plot
# retrieve data from DotPlot
plotData <- plot$data
# create two y axes
# define the clusters for the right y axis
plotData$id <- factor(plotData$id, levels = rev(c("Perineurial", "Chalice", "Fenestrated", "Subperineurial",
                                                  "Pseudo-cartridge", "Chiasm", "Cortex_1", "Cortex_2",
                                                  "Distal satellite", "Proximal satellite", "Epithelial",
                                                  "Ensheathing_1", "Ensheathing_2","Ensheathing_3",
                                                  "Marginal","Astrocyte_1","Astrocyte_2")))
ylabs2 <- levels(plotData$id)
# define the clusters for the left y axis
ylabs1 <- as.character(c(16,3,5,7,9,2,0,14,6,13,11,10,12,15,8,4,1))
# plot the data from DotPlot in graph with 2 y axis
final.plot <- plotData %>% filter(avg.exp.scaled > 0) %>% ggplot(aes(x = features.plot, y = as.numeric(id), color = avg.exp.scaled, size = pct.exp)) +
  geom_point() +
  scale_y_continuous( breaks = 1:length(ylabs1),
                      labels = ylabs1,
                      sec.axis = sec_axis(trans = ~.,
                                          breaks = 1:length(ylabs2),
                                          labels = ylabs2, name="Annotated Clusters",)) +
  scale_x_discrete() + theme_classic(14) +
  theme(axis.text.x = element_text(face = "italic", color="black",angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(size = 14, color = "black")) +
  labs(x= "Marker genes", y = "Clusters", color = "Avg. Exp", size = "% Exp") +
  scale_color_gradientn(colours=c("white", "#153743")) +
  scale_size_continuous(range=c(2, 6))
final.plot
ggsave('DotPlot_adult.png', final.plot, dpi = 900, device = "png")




#P15-----

exprMat <- read.csv("FlyphoneDB/ozel_p15/p15_matrix.csv", row.names = 1, check.names = FALSE)
cellInfo <- read.csv("FlyphoneDB/ozel_p15/p15_meta.csv", row.names = 1)
cellInfo$celltype <- as.character(cellInfo$celltype)
exprMat <- exprMat[ , row.names(cellInfo)]
# create seuratObj
seuratObj <- CreateSeuratObject(counts = exprMat)
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(seuratObj)
seuratObj <- ScaleData(seuratObj, features = all_genes)
seuratObj$celltype <- as.factor(cellInfo$celltype)
Idents(seuratObj) <- "celltype"
lamina_made<-subset(x = seuratObj, idents = c(107:110,117,257))




my_levels <- c(117,107,109,110,108,257)


lamina_made@active.ident <- factor(x = lamina_made@active.ident, levels = my_levels)

cluster_name <- c ("L5","L4","L1","L3","L2","LPC")
names(cluster_name) <- levels(lamina_made)

lamina_made<- RenameIdents(object = lamina_made, cluster_name)

df <- subset(Pathway_core_components, Pathway == 'HEDGEHOG signaling pathway')
genes <- df$Gene.Symbol

p15_plot <- DotPlot(lamina_made, scale= T,features = genes)+ theme(axis.text.x = element_text(angle = 90))
p15_plot


#Ramya's plot------
P15<-readRDS("~/FlyphoneDB/ozel_p15/GSE142787_P15.rds")

pbmc <- UpdateSeuratObject(object = P15)

rm(P15)

r_obj<-subset(x = pbmc, idents = c(91,253:256))
DefaultAssay(r_obj) <- "RNA"
r_gene<-c("numb",'dpn','elav')
DotPlot(r_obj, scale= T,features = r_gene)+ theme(axis.text.x = element_text(angle = 90))
#time plot----
## https://github.com/satijalab/seurat/issues/5117
library(SeuratData)
data("pbmc3k") 

# remove cells with NA in annotations and add a group column
pbmc3k <- subset(pbmc3k, cells = Cells(pbmc3k)[which(!is.na( pbmc3k$seurat_annotations))])
pbmc3k$group <- sample(x = c("t1", "t2" ), size = ncol(pbmc3k),  replace = T)
pbmc3k <- NormalizeData(pbmc3k)
Idents(pbmc3k) <- "seurat_annotations"
p1<- DotPlot(pbmc3k, features = "GNLY"  ) 

# make a new assay with the splited genes expression
GNLY.t1 <- pbmc3k[['RNA']]@data["GNLY",] * (pbmc3k$group == "t1")
GNLY.t2 <- pbmc3k[['RNA']]@data["GNLY",] * (pbmc3k$group == "t2")

pbmc3k[['NEW']] <- CreateAssayObject(data = rbind(GNLY.t1, GNLY.t2))
DefaultAssay(pbmc3k) <- "NEW"
p2 <- DotPlot(pbmc3k, features = c("GNLY.t1", "GNLY.t2")) 
p1+p2