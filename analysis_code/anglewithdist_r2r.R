rm(list = ls())#清空输出
PATH = getwd()
dir = paste(PATH,"/data/",sep="")
setwd(PATH)
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

matchtable_11 = c(0)
matchtable_12 = c(0)
matchtable_13 = c(0)
matchtable_22 = c(0)
matchtable_23 = c(0)
matchtable_33 = c(0)
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
  match_list_11 <- matrix(match_list, nrow = x)
  match_list_12 <- match_list_13 <- match_list_22 <- match_list_23 <- match_list_33 <- match_list_11
  
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
          match_list_11[j,k] = anglediff
        }
        
        if((flat_dot_dist1>28&flat_dot_dist1<65)&
           (flat_dot_dist2>28&flat_dot_dist2<65)){
          match_list_22[j,k] = anglediff
        }
        
        if(flat_dot_dist1>65&flat_dot_dist2>65){
          match_list_33[j,k] = anglediff
        }
        
        if((flat_dot_dist1<28&(flat_dot_dist2>28&flat_dot_dist2<65))|
           (flat_dot_dist2<28&(flat_dot_dist1>28&flat_dot_dist1<65))){
          match_list_12[j,k] = anglediff
        }
        if((flat_dot_dist1<28&flat_dot_dist2>65)|
           (flat_dot_dist2<28&flat_dot_dist1>65)){
          match_list_13[j,k] = anglediff
        }
        if((flat_dot_dist2>65&(flat_dot_dist1>28&flat_dot_dist1<65))|
           (flat_dot_dist1>65&(flat_dot_dist2>28&flat_dot_dist2<65))){
          match_list_23[j,k] = anglediff
        }
      }
    }
  }
  choose_11 = which(match_list_11!=0)
  match_list_11_all = match_list_11[choose_11]
  matchtable_11 = c(matchtable_11, match_list_11_all)
  
  choose_12 = which(match_list_12!=0)
  match_list_12_all = match_list_12[choose_12]
  matchtable_12 = c(matchtable_12, match_list_12_all)
  
  choose_13 = which(match_list_13!=0)
  match_list_13_all = match_list_13[choose_13]
  matchtable_13 = c(matchtable_13, match_list_13_all)
  
  choose_22 = which(match_list_22!=0)
  match_list_22_all = match_list_22[choose_22]
  matchtable_22 = c(matchtable_22, match_list_22_all)
  
  choose_23 = which(match_list_23!=0)
  match_list_23_all = match_list_23[choose_23]
  matchtable_23 = c(matchtable_23, match_list_23_all)
  
  choose_33 = which(match_list_33!=0)
  match_list_33_all = match_list_33[choose_33]
  matchtable_33 = c(matchtable_33, match_list_33_all)
  
}
matchtable_11 = matchtable_11[-1]
matchtable_12 = matchtable_12[-1]
matchtable_13 = matchtable_13[-1]
matchtable_22 = matchtable_22[-1]
matchtable_23 = matchtable_23[-1]
matchtable_33 = matchtable_33[-1]


A = matchtable_22
#A = matchtable_11
#A = matchtable_12
#A = matchtable_13
#A = matchtable_22
#A = matchtable_23
#A = matchtable_33

A[which(A>=90)] <- 180-A[which(A>=90)]

hist(A,breaks = seq(0,100,1),prob = TRUE, main = "parallel")

lines(density(A), col = "red")

















