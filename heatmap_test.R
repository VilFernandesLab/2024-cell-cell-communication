########
# Heat map test--------
exprMat <- read.csv("FlyphoneDB/online_sample/FlyPhone_Matrix_Example.csv", row.names = 1, check.names = FALSE)
cellInfo <- read.csv("~/FlyphoneDB/online_sample/FlyPhone_Metadata_Example.csv", row.names = 1)
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

avgexp <- AverageExpression(seuratObj, assay = "RNA", return.seurat = TRUE)


Pathway_core_components <- read.csv("FlyphoneDB/parameters/fly_phone_pathway_core_components.csv")

unique(Pathway_core_components$Pathway)
df <- subset(Pathway_core_components, Pathway == 'EGFR signaling pathway')
genes <- df$Gene.Symbol


DoHeatmap(avgexp, features = genes, label = TRUE ,draw.lines = FALSE, raster = FALSE, angle = 90) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))

# cluster re-assignment occurs, which re-assigns clustering in my_levels (assuming you have 12 clusters in total)
my_levels <- c("aEC1" ,   "aEC2"  ,  'aEC3' )

# Re-level object@ident
avgexp@active.ident <- factor(x = avgexp@active.ident, levels = my_levels)

#p15------
## actual seurat object----
data<- readRDS('~/FLyphoneDB/ozel_p15/GSE142787_P15.rds')

pbmc <- UpdateSeuratObject(object = data)

rm(data)

lamina<-subset(x = pbmc, idents = c(107:110,117,257))

Pathway_core_components <- read.csv("FlyphoneDB/parameters/fly_phone_pathway_core_components.csv")
unique(Pathway_core_components$Pathway)
df <- subset(Pathway_core_components, Pathway == 'NOTCH signaling pathway')
genes <- df$Gene.Symbol
avgexp <- AverageExpression(lamina, assay = "RNA", return.seurat = TRUE)



DoHeatmap(lamina, features = genes, label = TRUE ,draw.lines = TRUE, raster = FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))


## made seurat object----

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
avgexp <- AverageExpression(lamina_made, assay = "RNA", return.seurat = TRUE)



my_levels <- c(257,108,110,109,107,117)


lamina_made@active.ident <- factor(x = lamina_made@active.ident, levels = my_levels)

cluster_name <- c('LPC','L2','L3','L1','L4',"L5")
names(cluster_name) <- levels(lamina_made)

lamina_made<- RenameIdents(object = lamina_made, cluster_name)

DoHeatmap(lamina_made, features = genes, label = TRUE ,draw.lines = TRUE, raster = FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))

avgexp <- AverageExpression(lamina_made, assay = "RNA", return.seurat = TRUE)
DoHeatmap(avgexp, features = genes, label = FALSE ,draw.lines = FALSE, raster = FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))


#Adult------
exprMat <- read.csv("FlyphoneDB/ozel_adult/matrix.csv", row.names = 1, check.names = FALSE)
cellInfo <- read.csv("FlyphoneDB/ozel_adult/ad_meta.csv", row.names = 1)
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
lamina_made<-subset(x = seuratObj, idents = c(108,110,109,107,117))

my_levels <- c(108,110,109,107,117)

# Re-level object@ident
lamina_made@active.ident <- factor(x = lamina_made@active.ident, levels = my_levels)

cluster_name <- c('L2','L3','L1','L4',"L5")
names(cluster_name) <- levels(lamina_made)

lamina_made<- RenameIdents(object = lamina_made, cluster_name)

DoHeatmap(lamina_made, features = genes, label = TRUE ,draw.lines = TRUE, raster = FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))




avgexp <- AverageExpression(lamina_made, assay = "RNA", return.seurat = TRUE)
DoHeatmap(avgexp, features = genes, label = FALSE ,draw.lines = FALSE, raster = FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu")))



