start_time <- Sys.time()
set.seed(123)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(future.apply)
  library(Seurat)
  library(RColorBrewer)
  library(reshape2)
  library(network)
  library(igraph)
})

option_list = list(
  make_option(c("-i", "--matrix"), type="character", default=NULL,
              help="input matrix", metavar="character"),
  make_option(c("-a", "--metadata"), type="character", default=NULL,
              help="input metadata", metavar="character"),
  make_option(c("-p", "--lrpair"), type="character", default=NULL,
              help="annotation ligand receptor", metavar="character"),
  make_option(c("-s", "--subset"), type="character", default=NULL,
              help="input matrix", metavar="character"),
  make_option(c("-s", "--corecomponents"), type="character", default=NULL,
              help="annotation core components", metavar="character"),
  make_option(c("-c", "--cores"), type="character", default=NULL,
              help="number of cores", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output directory name", metavar="character"),
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(paste0("matrix input: ", opt$matrix))
print(paste0("metadata input: ", opt$metadata))
print(paste0("L-R pairs input: ", opt$lrpair))
print(paste0("ligend subset: ", opt$subset))
print(paste0("core components input: ", opt$corecomponents))
print(paste0("number of cores: ", opt$cores))
print(paste0("output directory: ", opt$output))

output_dir <- 'FlyphoneDB'
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# plan(multiprocess, workers = 8) ## => parallelize on your local computer
plan(multisession, workers = as.numeric(opt$cores)) ## => parallelize on your local computer

#install older version of future (>= 1.32.0)

####################################
# Input and cluster means
####################################


exprMat <- read.csv(opt$matrix, row.names = 1, check.names = FALSE)

cellInfo <- read.csv(opt$metadata, row.names = 1)
cellInfo$celltype <- as.character(cellInfo$celltype)


exprMat <- exprMat[ , row.names(cellInfo)]

exprMat <- sweep(exprMat, 2, Matrix::colSums(exprMat), FUN = "/") * 10000


LR_pairs <- read.csv(file = opt$lrpair, sep = ",")


gene_list <- unique(c(LR_pairs$Gene_secreted, LR_pairs$Gene_receptor))
common_genes <- intersect(gene_list, row.names(exprMat))

LR_pairs <- subset(LR_pairs, Gene_secreted %in% common_genes & Gene_receptor %in% common_genes)

exprMat <- as.matrix(exprMat)
exprMat <- t(exprMat)




df_Ligand <- exprMat[ , unique(LR_pairs$Gene_secreted)]
df_Receptor <- exprMat[ , unique(LR_pairs$Gene_receptor)]


celltype_df_Ligand <- cbind(cellInfo[ , c("celltype"), drop = FALSE], df_Ligand)
celltype_df_Receptor <- cbind(cellInfo[ , c("celltype"), drop = FALSE], df_Receptor)

df_group_by_celltype_Ligand <- celltype_df_Ligand %>%
  group_by(celltype) %>%
  summarise_all(mean) %>%
  as.data.frame()

row.names(df_group_by_celltype_Ligand) <- df_group_by_celltype_Ligand$celltype
df_group_by_celltype_Ligand$celltype <- NULL
df_group_by_celltype_Ligand <- t(df_group_by_celltype_Ligand)

# average Receptor counts by each celltype
df_group_by_celltype_Receptor <- celltype_df_Receptor %>%
  group_by(celltype) %>%
  summarise_all(mean) %>%
  as.data.frame()

row.names(df_group_by_celltype_Receptor) <- df_group_by_celltype_Receptor$celltype
df_group_by_celltype_Receptor$celltype <- NULL
df_group_by_celltype_Receptor <- t(df_group_by_celltype_Receptor)


####################################
# Interaction score-------
####################################

ligand_avg <- df_group_by_celltype_Ligand[LR_pairs$Gene_secreted, ] %>% as.data.frame()
receptor_avg <- df_group_by_celltype_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()

interaction_list <- list()
LR_pairs_one <- LR_pairs # combine

#for (i in sort(unique(cellInfo$celltype)) ) {
  LR_pairs_combine <- LR_pairs 
  #args<-commandArgs(TRUE)
  #i<-as.character(args[1]) # Array parameter
  
  i<-as.character(opt$subset)
  
  for (j in sort(unique(cellInfo$celltype)) ) {
    print(paste0(i, ">", j))
    LR_pairs_tmp <- LR_pairs
    LR_pairs_tmp[[paste0(i, ">", j, "_score")]] <- log1p(ligand_avg[[i]]) * log1p(receptor_avg[[j]])
    
    # permutation -------------------------------------------------------------
    # start permutatioin
    permutation_times <- 1000
    y <- future_lapply(1:permutation_times, function(ii) {
      cellInfo_sample_Ligand <- cellInfo
      cellInfo_sample_Ligand$celltype <- sample(cellInfo_sample_Ligand$celltype)
      
      celltype_df_sample_Ligand <- cbind(cellInfo_sample_Ligand[ , c("celltype"), drop = FALSE], df_Ligand)
      
      df_group_by_celltype_sample_Ligand <- celltype_df_sample_Ligand %>%
        group_by(celltype) %>%
        summarise_all(mean) %>%
        as.data.frame()
      
      row.names(df_group_by_celltype_sample_Ligand) <- df_group_by_celltype_sample_Ligand$celltype
      df_group_by_celltype_sample_Ligand$celltype <- NULL
      df_group_by_celltype_sample_Ligand <- t(df_group_by_celltype_sample_Ligand)
      
      # sample Receptor
      cellInfo_sample_Receptor <- cellInfo
      cellInfo_sample_Receptor$celltype <- sample(cellInfo_sample_Receptor$celltype)
      
      celltype_df_sample_Receptor <- cbind(cellInfo_sample_Receptor[ , c("celltype"), drop = FALSE], df_Receptor)
      
      df_group_by_celltype_sample_Receptor <- celltype_df_sample_Receptor %>%
        group_by(celltype) %>%
        summarise_all(mean) %>%
        as.data.frame()
      
      row.names(df_group_by_celltype_sample_Receptor) <- df_group_by_celltype_sample_Receptor$celltype
      df_group_by_celltype_sample_Receptor$celltype <- NULL
      df_group_by_celltype_sample_Receptor <- t(df_group_by_celltype_sample_Receptor)
      
      
      ####################################
      # Interaction score
      ####################################
      ligand_avg_tmp <- df_group_by_celltype_sample_Ligand[LR_pairs$Gene_secreted, ] %>% as.data.frame()
      receptor_avg_tmp <- df_group_by_celltype_sample_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()
      tmp <- log1p(ligand_avg_tmp[[i]]) * log1p(receptor_avg_tmp[[j]])
      tmp
    }, future.seed = TRUE)
    
    
    df <- data.frame(matrix(unlist(y), nrow=length(y), byrow=TRUE))
    df <- t(df)
    LR_pairs_tmp <- cbind(LR_pairs_tmp, df)
    
    
    LR_pairs_tmp$result <- rowSums(sapply(LR_pairs_tmp[, 13:ncol(LR_pairs_tmp)], function(x) x > LR_pairs_tmp[[paste0(i, ">", j, "_score")]]))
    LR_pairs_tmp[[paste0(i, ">", j, "_pvalues")]] <- LR_pairs_tmp$result / permutation_times
    
    
    LR_pairs_tmp <- LR_pairs_tmp[ , c(1:12, ncol(LR_pairs_tmp))]
    LR_pairs_tmp[LR_pairs_tmp[[paste0(i, ">", j, "_score")]] == 0, paste0(i, ">", j, "_pvalues")] <- 1
    
    LR_pairs_combine <- cbind(LR_pairs_combine,
                              LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
    )
    LR_pairs_one <- cbind(LR_pairs_one,
                          LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
    )
  }
  
  interaction_list[[i]] <- LR_pairs_combine
  


write.csv(LR_pairs_one, file = paste0(output_dir, "/", "interaction_list_",i,".csv"))

end_time <- Sys.time()

print(end_time - start_time)