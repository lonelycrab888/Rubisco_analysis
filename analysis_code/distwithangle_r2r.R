rm(list = ls())#清空输出
PATH = getwd()
dir_fun = paste(PATH,"/analysis_code/",sep="")
dir = paste(dir_fun,"/data/",sep="")
setwd(dir_fun)
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

matchtable = c(0)
matchtable_15 = c(0)
matchtable_45 = c(0)
matchtable_60 = c(0)
for (i in 1:numbermax){
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
  match_list <- rep(0,x*x)
  match_list <- matrix(match_list, nrow = x)
  match_list_15 <- match_list_45 <- match_list_60 <- match_list 
  
  for (j in 1:x){
    for(k in 1:x){
      if(E_mat_new[j,k]!=0){
        v_j = solveequ(M_new, j)
        v_k = solveequ(M_new, k)
        anglediff = as.numeric(findanglexy(v_j, v_k))
        flat_dot_dist1 = dist_flat_dot(M_new, j, k)
        flat_dot_dist2 = dist_flat_dot(M_new, k, j)
        deta = abs(flat_dot_dist1-flat_dot_dist2)
        # if ((flat_dot_dist1<28&flat_dot_dist2<28)|
        #    (flat_dot_dist2>65&flat_dot_dist1>65)){
        #   angle_list[j,k] = anglediff
        #  }
        match_list[j,k] = flat_dot_dist1
        if (anglediff<15|anglediff>165)
          match_list_15[j,k] = flat_dot_dist1
        if (anglediff<45|anglediff>135)
          match_list_45[j,k] = flat_dot_dist1
        if (anglediff<60|anglediff>120)
          match_list_60[j,k] = flat_dot_dist1
        
      }
    }
  }

  choose = which(match_list!=0)
  match_list_all = match_list[choose]
  matchtable = c(matchtable, match_list_all)
  
  choose_15 = which(match_list_15!=0)
  match_list_15_all = match_list_15[choose_15]
  matchtable_15 = c(matchtable_15, match_list_15_all)
  
  choose_45 = which(match_list_45!=0)
  match_list_45_all = match_list_45[choose_45]
  matchtable_45 = c(matchtable_45, match_list_45_all)
  
  choose_60 = which(match_list_60!=0)
  match_list_60_all = match_list_60[choose_60]
  matchtable_60 = c(matchtable_60, match_list_60_all)
  
}
matchtable = matchtable[-1]
matchtable_15 = matchtable_15[-1]
matchtable_45 = matchtable_45[-1]
matchtable_60 = matchtable_60[-1]

A = matchtable
A = matchtable_15
A = matchtable_45
A = matchtable_60


hist(A,breaks = seq(0,120,1),prob = TRUE, main = "<60")

lines(density(A), col = "red")

















