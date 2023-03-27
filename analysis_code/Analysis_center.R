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

##取较为完整的羧酶体进行整体空间的分析
#主要分析了rubisco在羧酶体中的分布规律以及朝向规律
Mnew <- M[which(M$TomoNumber==38),]
MM = findcenter(Mnew)
center = c(MM[1,2],MM[1,3],MM[1,4])
M_new = MM[-1,]

plotM = plotcentermap(center, M_new)
#plotM输出两列，
#第一列为空间朝向（实际z轴与拟合法向夹角）
#第二列为空间分布（实际rubisco与拟合羧酶体中心距离）
A = plotM[,1]
B = plotM[,2]
plot(plotM)
hist(A,breaks = seq(0,200,1),prob = TRUE, main = "parallel")
lines(density(A), col = "red")
hist(B,breaks = seq(0,1200,1),prob = TRUE, main = "distance2center")
lines(density(B),col = 'white')

plot(density(A))
plot(density(B))
