# general result process#

# Because flyphoneDB analyses were run in ARRAY JOB, this means the raw results were breaking down to
# about 200 files per stage (same number of the number of the clusters in that stage)
# The general workflow for this script:
# 1. Extract all putative interaction in each result file
# 2. Merge them into a big data frame 

# library -------
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(reshape2)
  library(dplyr)
})




#functions----------

#function to extract ligand cluster number 
call_l= function(x){
  return(str_split(x,'>')[[1]][1])
}

#function to extract receptor cluster number 
call_r= function(x){
  return(str_split(x,'>')[[1]][2])
}

#function to call cluster annotation
dict_df<-read.csv("FlyphoneDB/parameters/Clusters-dict.csv") #import dictionary file 
dict_df$Assignment <- gsub("*","",as.character(dict_df$Assignment),fixed = TRUE) #OPTIONAL: remove the * for the less confident assignment
call_cluster_name = function (x){
  name<-dict_df[which(dict_df$Cluster == x),2]
  return(name)
}

#function to extract putative results from each result file
make_clean_df = function (df){
  clu_sqr=(ncol(df)-3)/2
  clean_df<-data.frame(matrix(ncol = 5, nrow = 0))
  colnames(clean_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")
  for (x in c(1:clu_sqr)) {
    s=(x-1)*2+4
    p=(x-1)*2+5
    
    fil_df<- df%>%
      filter(
        .[[p]] <= 0.5 ##DEFINE P-VALUE HERE!!!!!!!!! ----
      ) 
    
    if (nrow(fil_df)==0){
      
    }else{
      fil_df<-fil_df[,c(1:3,s)]
      colnames(fil_df)[ncol(fil_df)]<-'score'
      fil_df['pairing']<-str_split(colnames(df)[s],'_')[[1]][1]
      
      clean_df<-rbind(clean_df,fil_df)
    }
  }
  
  clean_df$l_cluster<-lapply(clean_df$pairing,call_l)
  clean_df$l_cluster<-lapply(clean_df$l_cluster, call_cluster_name)
  
  clean_df$r_cluster<-lapply(clean_df$pairing,call_r)
  clean_df$r_cluster<-lapply(clean_df$r_cluster, call_cluster_name)
  
  return(clean_df)
}

# read in------

## CHANGE STAGE HERE!!!!!!!!!!----
stage='p15' 


meta <- read.csv(paste0("FlyphoneDB/ozel_",stage,"/meta.csv")) # input meta data 
celltype<-unique(meta$celltype) # extract cell type
file<-list.files(path = paste0("FlyphoneDB/ozel_",stage,"/output_test_dataset")) # save all file path
length(file) # read out of how many output files you have 

#Actual process----

#build an empty data frame to save all your output 
big_df<-data.frame(matrix(ncol = 5, nrow = 0))
colnames(big_df)<-c("Gene_secreted","Gene_receptor","pathway_receptor","score","pairing")

#iterate through all the files
for (f in file){
  path<-paste0("FlyphoneDB/ozel_",stage,"/output_test_dataset/",f)
  df_o<-read.csv(path,check.names = FALSE)
  df<- data.frame(df_o[,-1], row.names = df_o[,1],check.names = FALSE)
  c_df<-make_clean_df(df) # extract putative results
  big_df<-rbind(big_df,c_df) # append to the big data frame 
  print(f) #print out path name just to check progress of iteration 
}

#filter out interaction with cluster 0 (junk) 
big_df<-big_df %>%
  filter(r_cluster!="character(0)")
big_df<-big_df %>%
  filter(l_cluster!="character(0)")

#formatting 
big_df$l_cluster <- sapply(big_df$l_cluster, as.character)
big_df$r_cluster <- sapply(big_df$r_cluster, as.character)

#write out 
write.csv(big_df, file = paste0("FlyphoneDB/ozel_",stage,"/all_p0.05.csv"),row.names = FALSE)



#read in (for entry length)------
stage<-'p15'
df<-read.csv(paste0("FlyphoneDB/ozel_",stage,"/all_p0.05.csv"))
nrow(df)
