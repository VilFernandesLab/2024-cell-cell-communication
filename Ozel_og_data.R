#######################
### Ozel et al. data ##
#######################

library(R.utils)
library(Seurat)

gunzip("GSE142787_P15.rds.gz", remove=FALSE)

data<- readRDS('~/FLyphoneDB/ozel_adult/GSE142787_Adult.rds')

pbmc <- UpdateSeuratObject(object = data)

rm(data)

table(pbmc$orig.ident)


# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


pbmc <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)



# Basic Features-------

length(Cells(pbmc)) # number of cells: 109743
length(rownames(pbmc)) # number of genes: 2000
Layers(pbmc)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)


pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)


pbmc <- RunUMAP(pbmc, dims = 1:10)


DimPlot(pbmc, reduction = "tsne", group.by = "FinalIdents")


cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# subset------

## L3 subset------

l3<-subset(x = pbmc, subset = L3_vs_P15 =='2')

# Randomly selected 3 (including 1 glia cluster) and try to convert to flyphone DB format for web

g_example<-subset(x = pbmc, idents = c(198,1,3))

saveRDS(g_example, file = "~/FlyphoneDB/1,3,198_subset/1,3,198.rds")

write.csv(g_example@active.ident, file='~/FlyphoneDB/gex_meta.csv', quote=FALSE, sep='\t', col.names = TRUE)
write.csv(g_example@assays[["RNA"]]@counts, file='sample_Gene_Count_per_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)

## meta---
gex_meta<-data.frame(g_example@active.ident)
colnames(gex_meta)<-'celltype'
write.csv(gex_meta, file='~/FlyphoneDB/gex_meta.csv', 
          quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)

gex_meta<- read.csv("~/FlyphoneDB/gex_meta.csv")

##matrix-----

gex_matrix<-data.frame(g_example@assays[["RNA"]]@counts)
write.csv(gex_matrix, file='~/FlyphoneDB/gex_matrix.csv', 
          quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)

gex_matrix <- read.csv("~/FlyphoneDB/1,3,198_subset/gex_matrix.csv")

DimPlot(g_example, reduction = "tsne", group.by = "FinalIdents")
FeaturePlot(g_example, features = c("dlp",'dally')) # high 1>198 value
FeaturePlot(g_example, features = 'Dscam2')


#Exporting Function----

export_m<- function (seuobj,meta_path,matrix_path){
  meta<-data.frame(seuobj@active.ident)
  colnames(meta)<-'celltype'
  write.csv(meta, file=meta_path, 
            quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)
  
  matrix<-data.frame(seuobj@assays[["RNA"]]@counts)
  write.csv(matrix, file=matrix_path, 
            quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)
}

export_m(pbmc,"~/FlyphoneDB/ozel_adult/ad_meta.csv","~/FlyphoneDB/ozel_adult/ad_matrix.csv")


p15_matrix<- read.csv("~/FlyphoneDB/ozel_adult/ad_matrix.csv")




# UMAP-----
DefaultAssay(pbmc)<-'RNA'

pbmc_small <- RunUMAP(object = l3, dims = 1:150)
FeaturePlot(pbmc, features = c('ple'),label=TRUE)
