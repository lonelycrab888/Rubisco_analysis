rm(list = ls())#清空输出
PATH = getwd()
setwd(PATH)
dir = paste(PATH,"/data/",sep="")

maxdist = 110


library(xlsx)
library(readxl)
library(tidyverse)
library(factoextra)
library(cluster)
library(rgl)
library("MASS")
library("ggpubr")
library(ggplot2)
library(fpc)
source('function.R')

M <- read.table(paste(dir,"matching_checknew.xls",sep=""),head=FALSE)
colnames(M) <- c("MicrographName","CoordinateX","CoordinateY",
                 "CoordinateZ","AngleRot","AngleTilt","Angle",
                 "TomoNumber")
M$CoordinateX <- round(M$CoordinateX, 2) 
M$CoordinateY <- round(M$CoordinateY, 2)
M$CoordinateZ <- round(M$CoordinateZ, 2)
M$Angle <- round(M$Angle,3)
numbermax <- M[nrow(M),8]

matchtable_hth = c(0)
matchtable_hts = c(0)
matchtable_sts = c(0)
matchtable_other = c(0)

for (i in 1:numbermax){
  i = 1
  Mnew <- M[which(M$TomoNumber==i),]
  distxy <- Mnew[,2:3]
  M_cluster <- clustering_xy(distxy , Mnew)
  
  #l_index:聚类个数
  l_index = M_cluster[1,9]
  
  #flag_cluster:聚类组分，当组分只有0时值为0，否则为-1
  flag_cluster = M_cluster[1,10]
  
  #判断是否聚类成球
  if(l_index==0 | (l_index==1 & flag_cluster==0)){
    next
  }
  
  #M_cluster:聚类后，rubisco集合
  M_cluster <- M_cluster[,1:8]
  distancedata_cluster <- M_cluster[,2:4]
  #构建距离矩阵E_mat
  E_mat <- as.matrix(dist(distancedata_cluster,p=2))
  #将两两距离大于maxdist的元素置0
  E_mat[E_mat > maxdist] <- 0
  useless_dis <- which(rowSums(E_mat)==0)
  #M_new:有互作的rubisco集合
  M_new <- M_cluster[-useless_dis,]
  distancedata_new <- M_new[,2:4]
  
  #构建新的有互作的rubisco的距离矩阵
  E_mat_new <- as.matrix(dist(distancedata_new,p=2))
  E_mat_new[E_mat_new>maxdist] <- 0
  
  x = nrow(E_mat_new)
  
  #构建全0空矩阵，记录有互作的rubisco间感兴趣的参数
  match <- rep(0,x*x)
  match <- matrix(match_list, nrow = x)
  Match_hts <- Match_hth <- Match_sts <- Match_other <- match
  
  
  for (j in 1:x){
    for(k in j:x){
      if(E_mat_new[j,k]!=0){
        v_j = solveequ(M_new, j)
        v_k = solveequ(M_new, k)
        anglediff = as.numeric(findanglexy(v_j, v_k))
        flat_dot_dist1 = dist_flat_dot(M_new, j, k)
        flat_dot_dist2 = dist_flat_dot(M_new, k, j)
        deta = abs(flat_dot_dist1-flat_dot_dist2)
        if(flat_dot_dist1<28&flat_dot_dist2<28){
          if(anglediff<35|anglediff>145){
            Match_sts[j,k] = 1
          }
          else Match_other[j,k] = 1
        }
        
        else if(flat_dot_dist1>65&flat_dot_dist2>65){
          if(anglediff<35|anglediff>145){
            Match_hth[j,k] = 1
          }
          else Match_other[j,k] = 1
        }
        
        else if((flat_dot_dist1<28&flat_dot_dist2>65)|
           (flat_dot_dist2<28&flat_dot_dist1>65)){
          if(anglediff<115|anglediff>65){
            Match_hts[j,k] = 1
          }
          else Match_other[j,k] = 1
        }
        else Match_other[j,k] = 1
        
      }
    }
  }
  Match_hth = t(Match_hth)+Match_hth
  
  Match_hts = t(Match_hts)+Match_hts
  
  Match_sts = t(Match_sts)+Match_sts
  
  Match_other = t(Match_other)+Match_other
  
  M_hth= M_new[-which(rowSums(Match_hth)==0),]
  
  M_hts= M_new[-which(rowSums(Match_hts)==0),]
  
  M_sts = M_new[-which(rowSums(Match_sts)==0),]
  
  M_other = M_new[-which(rowSums(Match_other)==0),]
  
  matchtable_hth = rbind(matchtable_hth,M_hth)
  
  matchtable_hts = rbind(matchtable_hts,M_hts)
  
  matchtable_sts = rbind(matchtable_sts,M_sts)
  
  matchtable_other = rbind(matchtable_otehr,M_other)
  
}
matchtable_hth = matchtable_hth[-1]
matchtable_hts = matchtable_hts[-1]
matchtable_sts = matchtable_sts[-1]
matchtable_other = matchtable_other[-1]

write.xlsx(matchtable_hth, paste(dir,"hth_all.xls",sep=""), append = TRUE)
write.xlsx(matchtable_hts, paste(dir,"hts_all.xls",sep=""), append = TRUE)
write.xlsx(matchtable_sts, paste(dir,"sts_all.xls",sep=""), append = TRUE)
write.xlsx(matchtable_other, paste(dir,"other_all.xls",sep=""), append = TRUE)





















