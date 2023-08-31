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

#M<- read.table("F:/git-load/Rubisco_analysis/analysis_code/data/matching_checknew.xls",header = FALSE)

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
#A = matchtable_45
#A = matchtable_60


hist(A,breaks = seq(0,120,1),prob = TRUE, main = "<60")

lines(density(A), col = "red")


plot(density(A), # x和y坐标
     type = "l", # 图的类型
     main = "", # 图的主标题
     xlab = "Distance(Å)", ylab = "Density", # x、y轴标注
     ann = TRUE, # 逻辑值，是否使用默认的x、y轴标注注释
     axes = TRUE, # 逻辑值，是否显示坐标轴， "xaxt" 或 "yaxt" 选择不显示对应坐标轴
     #frame.plot = axes, # 是否显示图边框
     #asp = 10000, # y/x 的比例
     xgap.axis = 0.1, # x轴标签显示的距离
     ygap.axis = 0.1,# y轴标签显示的距离
     bty = 'o' # 图边框类型
     
)



A <- data.frame(sort(A))
colnames(A) <- "distance"



#pdf("C:/Users/zhr/Desktop/test.pdf",width=10,height=6)

P_density = ggplot(A, aes(x=distance))+ geom_density(size = 0.8)+ theme_classic()+ 
  scale_x_continuous(expand = c(0, 0),name = "Distance(Å)", limits = c(0,120),
                     breaks = c(0,30,60,90,120)) + 
  scale_y_continuous(expand = c(0, 0),name = "Probability Density", limits = c(0,0.02), 
                     breaks = c(0.005, 0.01,  0.015, 0.02))+
  theme(
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    plot.margin=unit(rep(1,4),'cm'))
print(P_density)
#dev.off()







length(density(A$distance)$x)

cdf_data <- data.frame(
  distance = density(A$distance)$x,
  cdf = cumsum(density(A$distance)$y * diff(density(A$distance)$x))
)

# Create plot for cumulative distribution function (CDF)
p_cdf <- ggplot(cdf_data, aes(x = distance, y = cdf)) +
  geom_line(size = 0.8) +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0), name = "Distance(Å)", limits = c(0, 120), breaks = c(0, 30, 60, 90, 120)) +
  scale_y_continuous(expand = c(0, 0), name = "Cumulative Probability", limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme(
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    plot.margin = unit(rep(1, 4), 'cm')
  )








