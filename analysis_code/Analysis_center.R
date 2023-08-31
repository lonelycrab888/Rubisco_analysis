rm(list = ls())#清空输出
PATH = getwd()
dir_fun = paste(PATH,"/analysis_code/",sep="")
dir = paste(dir_fun,"data/",sep="")
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
hist(A,breaks = seq(0,200,5),prob = TRUE, main = "parallel")
lines(density(A), col = "red")
hist(B,breaks = seq(0,1200,5),prob = TRUE, main = "distance2center")
lines(density(B),col = 'red')

plot(density(A), # x和y坐标
     type = "l", # 图的类型
     main = "", # 图的主标题
     xlab = "Angle(°)", ylab = "Density", # x、y轴标注
     ann = TRUE, # 逻辑值，是否使用默认的x、y轴标注注释
     axes = TRUE, # 逻辑值，是否显示坐标轴， "xaxt" 或 "yaxt" 选择不显示对应坐标轴
     #frame.plot = axes, # 是否显示图边框
     #asp = 10000, # y/x 的比例
     xgap.axis = 0.1, # x轴标签显示的距离
     ygap.axis = 0.1,# y轴标签显示的距离
     bty = 'o' # 图边框类型
     
)

plot(density(B), # x和y坐标
     type = "l", # 图的类型
     main = "", # 图的主标题
     xlab = "Distance(Å)", ylab = "Density", # x、y轴标注
     ann = TRUE, # 逻辑值，是否使用默认的x、y轴标注注释
     axes = TRUE, # 逻辑值，是否显示坐标轴， "xaxt" 或 "yaxt" 选择不显示对应坐标轴
     #frame.plot = axes, # 是否显示图边框
     #asp = 1, # y/x 的比例
     xgap.axis = 0.1, # x轴标签显示的距离
     ygap.axis = 0.1,# y轴标签显示的距离
     bty = 'L' # 图边框类型
     
)

B <- data.frame(sort(B))
colnames(B) <- "distance"



pdf("C:/Users/zhr/Desktop/test2.pdf",width=4,height=4)

ggplot(B, aes(x=distance))+ geom_density(size = 0.8)+ theme_classic2()+ 
  scale_x_continuous(expand = c(0,0),name = "Distance(Å)", limits = c(0,1400), breaks = c(0,200,400,600,800,1000,1200,1400)) + 
  scale_y_continuous(expand = c(0,0),name = "Density", limits = c(0,0.0028), breaks = c(0.0005, 0.001, 0.0015, 0.002, 0.0025))+
  theme(
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    plot.margin=unit(rep(1,4),'cm'))


dev.off()




A <- data.frame(sort(A))
colnames(A) <- "Angle"

pdf("C:/Users/zhr/Desktop/test.pdf",width=4,height=4)

ggplot(A, aes(x=Angle))+ geom_density(size = 0.8)+ theme_classic2()+ 
  scale_x_continuous(expand = c(0,0),name = "Angle(°)",limits = c(0,200), breaks = c(0,45, 90,135,180)) + 
  scale_y_continuous(expand = c(-0.002,-0.002),name = "Density",limits = c(0,0.0095), breaks = c(0.002,0.003,0.004, 0.005,0.006,0.007))+
  theme(
    axis.title.x = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    plot.margin=unit(rep(1,4),'cm'))

dev.off()


