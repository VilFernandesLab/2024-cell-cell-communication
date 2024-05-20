## Make matrix and metadata##

library(Seurat)
library(R.utils)

#input-----



#gunzip("FlyphoneDB/ozel_l3/GSE167266_OL.L3_P15_merged.rds.gz", remove=FALSE)

data<- readRDS('FlyphoneDB/ozel_p50/GSE142787_P50.rds')

pbmc <- UpdateSeuratObject(object = data)

rm(data)

pbmc <- RenameCells(object = pbmc, add.cell.id = "X")


## meta-----
meta<-data.frame(pbmc@active.ident)
colnames(meta)<-'celltype'
write.csv(meta, file='FlyphoneDB/ozel_p50/meta.csv', 
          quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)



##matrix-----

matrix<-data.frame(pbmc@assays[["RNA"]]@counts)
write.csv(matrix, file='FlyphoneDB/ozel_p50/matrix.csv', 
          quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)




celltype<-unique(meta$celltype)
length(celltype)
print(celltype)
