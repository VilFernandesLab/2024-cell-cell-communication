library(R.utils)
library(Seurat)

#gunzip("GSE142787_P15.rds.gz", remove=FALSE)

#input----
data<- readRDS('~/FLyphoneDB/ozel_p15/GSE142787_p15.rds')

pbmc <- UpdateSeuratObject(object = data)

rm(data)

#rename cells----

pbmc <- RenameCells(object = pbmc, add.cell.id = "X")
head(x = colnames(x = pbmc))


# output function-----
export_m<- function (seuobj,meta_path,matrix_path){
  meta<-data.frame(seuobj@active.ident)
  colnames(meta)<-'celltype'
  write.csv(meta, file=meta_path, 
            quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)
  
  matrix<-data.frame(seuobj@assays[["RNA"]]@counts)
  write.csv(matrix, file=matrix_path, 
            quote=FALSE, sep='\t', row.names= TRUE,col.names = TRUE)
}

export_m(pbmc,"~/FlyphoneDB/ozel_p15/p15_meta.csv","~/FlyphoneDB/ozel_p15/p15_matrix.csv")
