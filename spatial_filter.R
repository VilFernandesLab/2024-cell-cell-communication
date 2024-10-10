# SPATIAL FILTER FOR ADULT STAGE #
# Based on Drosophila optic lobe cell types projections


library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)


# read in spatial parameter files

spatial<-read.csv("FlyphoneDB/parameters/spatial_out.csv")
spatial[is.na(spatial)] <- 0
spatial<- data.frame(spatial[,-1], row.names = spatial[,1])

#test-----

# test non overlap 'R1 and Mt5'
# test overlap R1 and R2


t_df<-spatial[c('L1',"Pm1"),]
t_list<-as.numeric(colSums(t_df))
o_df<-spatial[c('R1',"R2"),]
o_list<-as.numeric(colSums(o_df))

## overlap ----
if(2 %in% t_list == TRUE){
  print('overlap')
}else{
  print('non_overlap')
}


## section proximity-----

#Assumption:
# We only interested in the nueropil which has multiple layer ME,LO,LOP

p_df<-read.csv('FlyphoneDB/parameters/test_distant.csv')
p_df<- data.frame(p_df[,-1], row.names = p_df[,1])

# 1. test whether both cell types present in the same neuropil

#ME
1 %in% p_df[1,c(5:14)]&1 %in% p_df[2,c(5:14)]
#LO
1 %in% p_df[1,c(16:21)]&1 %in% p_df[2,c(16:21)]
#LOP
1 %in% p_df[1,c(22:25)]&1 %in% p_df[2,c(22:25)]
#FALSE: no close overlap

# 2. if both present (step 1 true), whether there is no layer gap between two cell type
sum_list<-as.numeric(colSums(p_df))
#ME
0 %in% sum_list[c(5:14)]
#LO
0 %in% sum_list[c(16:21)]
#LOP
0 %in% sum_list[c(22:25)]
#FLASE: no close overlap

# 3. if there is gap (step 2 true), where is the gap 

# 4. if the gap is in the middle of the neuropil layers: 
#    e.g. ME 3-4-5, then there is no overlapping because both cell types comes from either side of the neuropil
#    no overlap
#    but if the gap is at the most outside layers ( e.g. ME including layer 1 or 10)
#    There is a chance that one of the cell types situated in the middle of the neuropil like mt8
#    RESULT: rint out layer number 

# 5. (potential step)
# result with 7-8-9-10, test whether number are continuous?


# Functions------

# Function for overall overlapping for two cell types 
# TO filter out non overlapping 

## basic functions------
overall_overlap<-function(cell_a,cell_b){
  
  if(cell_a%in% rownames(spatial)&cell_b%in% rownames(spatial)){
    df<-spatial[c(cell_a,cell_b),]
    sm<-as.numeric(colSums(df))
    if(2 %in% sm == TRUE){
      return('overlap')
    }else{
      return('non_overlap')
    }
  }else{
    return('unkown')
  }
  
}




#below is test and proof of concept, no need to run 
brk<-c(1,1,1,0,0,0,1,1,0,0)
con<-c(1,1,1,1,1,1,0,0,0,0)

b_list<-which(brk==0)
b_test<-c(b_list[1]:b_list[length(b_list)])
length(b_list)==length(b_test)

c_list<-which(con==0)
c_test<-c(c_list[1]:c_list[length(c_list)])
length(c_list)==length(c_test)
# test end HERE


# function testing whether there is any neurophil layer gap between the two cell types
test_contin<- function (sm){
  c_list<-which(sm==0)
  c_test<-c(c_list[1]:c_list[length(c_list)])
  return(length(c_list)==length(c_test))
}
test_contin(sm[c(5:14)])


##neuropil functions-----
# infer spatial relationship in EACH neurophil 

# MEDULA
me_proximity<-function(cell_a,cell_b){
  df<-spatial[c(cell_a,cell_b),]
  sm<-as.numeric(colSums(df))
  if(1 %in% df[1,c(5:14)]==FALSE | 1 %in% df[2,c(5:14)]==FALSE){
    return('non_overlap')
  }else{
    if(0 %in% sm[c(5:14)]==FALSE){
      return('adjacent')
    }else{
      if(1 %in% which(sm[c(5:14)]==0)==TRUE | 10 %in% which(sm[c(5:14)]==0)==TRUE ) {  # if 0 located at the two outside layers of 
        if(test_contin(sm[c(5:14)])==TRUE){
          return('adjacent')
        }else{
          return('non_overlap')
        }
      }else{
        return('non_overlap')
      }
    }
  }
}


# LOBULA
lo_proximity<-function(cell_a,cell_b){
  df<-spatial[c(cell_a,cell_b),]
  sm<-as.numeric(colSums(df))
  if(1 %in% df[1,c(16:21)]==FALSE | 1 %in% df[2,c(16:21)]==FALSE){
    return('non_overlap')
  }else{
    if(0 %in% sm[c(16:21)]==FALSE){
      return('adjacent')
    }else{
      if(1 %in% which(sm[c(16:21)]==0)==TRUE | 6 %in% which(sm[c(16:21)]==0)==TRUE ) {
        if(test_contin(sm[c(16:21)])==TRUE){
          return('adjacent')
        }else{
          return('non_overlap')
        }
      }else{
        return('non_overlap')
      }
    }
  }
}


# LOBULA PLATE
lop_proximity<-function(cell_a,cell_b){
  df<-spatial[c(cell_a,cell_b),]
  sm<-as.numeric(colSums(df))
  if(1 %in% df[1,c(22:25)]==FALSE | 1 %in% df[2,c(22:25)]==FALSE){
    return('non_overlap')
  }else{
    if(0 %in% sm[c(22:25)]==FALSE){
      return('adjacent')
    }else{
      if(1 %in% which(sm[c(22:25)]==0)==TRUE | 4 %in% which(sm[c(22:25)]==0)==TRUE ) {
        if(test_contin(sm[c(22:25)])==TRUE){
          return('adjacent')
        }else{
          return('non_overlap')
        }
      }else{
        return('non_overlap')
      }
    }
  }
}





# try on------

score<-read.csv('FlyphoneDB/ozel_adult/all_p0.05.csv') 

score$overlap<-mapply(overall_overlap,score$l_cluster,score$r_cluster)


for (i in c(1:nrow(score))){
  if(score[i,8]=='non_overlap'){
    score$ME[i]<-me_proximity(score$l_cluster[i],score$r_cluster[i])
    score$LO[i]<-lo_proximity(score$l_cluster[i],score$r_cluster[i])
    score$LOP[i]<-lop_proximity(score$l_cluster[i],score$r_cluster[i])
  }else{
    score$ME[i]<-'N_A'
    score$LO[i]<-'N_A'
    score$LOP[i]<-'N_A'
  }
}

overlap<-score %>% filter(overlap=='overlap'|ME=='adjacent'|LO=='adjacent'|LOP=='adjacent') # including overlaps or adjacent
write.csv(overlap,file = "FlyphoneDB/ozel_adult/p0.05_spatial_fil.csv",row.names = FALSE)
# ME=='adjacent'|LO=='adjacent'|LOP=='adjacent'



