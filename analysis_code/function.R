
#构建旋转矩阵 欧拉角为（a,b,c）
EAZYZtoRM <- function(a,b,c){
  R11 = cos(a)*cos(b)*cos(c)-sin(a)*sin(c)
  R12 = cos(a)*sin(c)+cos(b)*cos(c)*sin(a)
  R13 = -cos(c)*sin(b)
  R21 = -cos(c)*sin(a)-cos(a)*cos(b)*sin(c)
  R22 = cos(a)*cos(c)-cos(b)*sin(a)*sin(c)
  R23 = sin(c)*sin(b)
  R31 = cos(a)*sin(b)
  R32 = sin(b)*sin(a)
  R33 = cos(b)
  RM = matrix(c(R11,R12,R13,R21,R22,R23,R31,R32,R33),ncol=3)
  return(RM)
}


#求索引为i的rubisco的中心坐标系的xy平面在基础坐标系中的平面方程:Ax+By+Cz+const = 0
solveequ <- function(M, i){
  x = M[i,2]
  y = M[i,3]
  z = M[i,4]
  a = M[i,5]
  b = M[i,6]
  c = M[i,7]
  a = a*pi/180
  b = b*pi/180
  c = c*pi/180
  
  RM = EAZYZtoRM(a,b,c)
  v = as.matrix(RM[3,],ncol = 1)
  v_xyz = matrix(c(x,y,z),nrow = 1)
  D = - v_xyz %*% v
  v_new = matrix(c(v[1,1],v[2,1],v[3,1], D),nrow=4)
  return(v_new)
}

#修正方程，使两平面平行
solveequ_correct <- function(M,i,j,vi,vj){
  v_correct = vj+vi
  v_correct_i = v_correct_j = v_correct
  xi = M[i,2]
  yi = M[i,3]
  zi = M[i,4]
  xj = M[j,2]
  yj = M[j,3]
  zj = M[j,4]
  vi_xyz = matrix(c(xi,yi,zi),nrow = 1)
  vj_xyz = matrix(c(xj,yj,zj),nrow = 1)
  vi_new = as.matrix(v_correct[-4,],ncol = 1)
  vj_new = as.matrix(v_correct[-4,],ncol = 1)
  Di = -vi_xyz %*% vi_new
  Dj = -vj_xyz %*% vj_new
  v_correct_i[4,1] = as.numeric(Di)
  v_correct_j[4,1] = as.numeric(Dj)
  v_new = cbind(v_correct_i,v_correct_j)
  return(v_new)
}


dist_flat <- function(v){
  p = abs(v[4,1]-v[4,2])
  q = sqrt(v[1,1]*v[1,1]+v[2,1]*v[2,1]+v[3,1]*v[3,1])
  A = p/q
  return(A)
}

dist_flat_dot <- function(M,i,j){
  v = solveequ(M, i)
  x = M[j,2]
  y = M[j,3]
  z = M[j,4]
  O = matrix(c(x,y,z,1),nrow = 1)
  p = O %*% v
  q = sqrt(v[1,1]*v[1,1]+v[2,1]*v[2,1]+v[3,1]*v[3,1])
  return(abs(p/q))
}


dist_line <- function(M_new,j,k,correct){
  x1 = M_new[j,2]
  y1 = M_new[j,3]
  z1 = M_new[j,4]
  x2 = M_new[k,2]
  y2 = M_new[k,3]
  z2 = M_new[k,4]  
  M1M2 = v1 = matrix(c(x1-x2,y1-y2,z1-z2),ncol = 1)
  v2 = matrix(c(correct[1,1],correct[2,1],correct[3,1]),ncol = 1)
  
  v = c(v1[2,1]*v2[3,1]-v1[3,1]*v2[2,1],
        v1[3,1]*v2[1,1]-v1[1,1]*v2[3,1],
        v1[1,1]*v2[2,1]-v1[2,1]*v2[1,1]
  )
  v = as.matrix(v,nrow = 1)
  p = sqrt(sum(v**2))
  q = sqrt(sum(v2**2))
  
  #p = crossprod(v,M1M2)
  #q = sqrt(sum(v**2))
  result = as.numeric(abs(p/q))
  return(result)
}

#求索引为i的rubisco的中心坐标系的xy平面在基础坐标系中的平面方程:Ax+By+Cz+const = 0
findanglexy <- function(v1, v2){
  v1 = as.matrix(v1[-4,],ncol = 1)
  v2 = as.matrix(v2[-4,],ncol = 1)
  dot = t(v1) %*% v2
  theta = dot / (sqrt(sum(v1**2))*sqrt(sum(v2**2)))
  theta = acos(theta)*180/pi
  return(theta)
}


#创建期待向量
create_expectV <- function(p,M,i){
  x <- M[i,2]
  y <- M[i,3]
  z <- M[i,4]
  xnew = x-p[1]
  ynew = y-p[2]
  znew = z-p[3]
  v_e = matrix(c(xnew,ynew,znew,0),ncol = 1)
  return(v_e)
}

center2point <- function(center, M, i){
  x <- M[i,2]
  y <- M[i,3]
  z <- M[i,4]
  x_c = center[1]
  y_c = center[2]
  z_c = center[3]
  dist = sqrt(((x-x_c)**2)+
                ((y-y_c)**2)+
                ((z-z_c)**2))
  return(dist)
}


clustering_xy <- function(distxy, Mnew){
  #聚类DBSCAN算法（对xy面投影进行）
  db = dbscan(distxy, eps = 128, MinPts = 5, scale = FALSE)
  #用于作聚类图
  #fviz_cluster(db,distxy,stand = FALSE, frame = FALSE, geom = "point")
  
  cluster <- db$cluster
  total <- c(sum(cluster==1),sum(cluster==2),sum(cluster==3),sum(cluster==4),
             sum(cluster==5),sum(cluster==6),sum(cluster==7),sum(cluster==8))
  index_total <- which(total>50)
  l_index <- length(index_total)
  cluster_logi <- cluster
  for(i in 1:l_index){
    cluster_logi[which(cluster==index_total[i])] <- 1
  }
  cluster_logi[which(cluster_logi!=1)] <- 0
  useless_cluster <- which(cluster_logi==0)
  M_cluster <- Mnew[-useless_cluster,]
  LL = -1
  if(l_index==1) {
    LL=index_total[1]
  }
  M_cluster_new <- cbind(M_cluster, l_index, LL)
  return(M_cluster_new)
}

clustering_xz <- function(distxz, M_cluster){
  dbxz = dbscan(distxz, eps = 120, MinPts = 5, scale = FALSE)
  
  clusterxz <- dbxz$cluster
  total <- c(sum(clusterxz==1),sum(clusterxz==2),sum(clusterxz==3),sum(clusterxz==4),
             sum(clusterxz==5),sum(clusterxz==6),sum(clusterxz==7),sum(clusterxz==8))
  index_total <- which(total>50)
  l_total <- length(index_total)
  cluster_logi <- clusterxz
  for(i in 1:l_total){
    cluster_logi[which(clusterxz==index_total[i])] <- 1
  }
  cluster_logi[which(cluster_logi!=1)] <- 0
  useless_cluster <- which(cluster_logi==0)
  M_cluster_new <- M_cluster[-useless_cluster,]
  return(M_cluster_new)
  
}


findcenter <- function(Mnew){
  distxy <- Mnew[,2:3]
  
  M_cluster <- clustering_xy(distxy , Mnew)
  l_index = M_cluster[1,9]
  flag_cluster = M_cluster[1,10]
  M_cluster <- M_cluster[,1:8]
  distancedata_cluster <- M_cluster[,2:4]
  y_bar = 0.5*(min(distancedata_cluster$CoordinateY)+
                 max(distancedata_cluster$CoordinateY))
  x_bar = 0.5*(min(distancedata_cluster$CoordinateX)+
                 max(distancedata_cluster$CoordinateX))
  distxz <- M_cluster[,2:4]
  distxz <- distxz[,-2]
  M_cluster_xyz <- clustering_xz(distxz , M_cluster)
  M_new = M_cluster_xyz
  distancedata_cluster_new <- M_new[,2:4]
  z_bar = 0.5*(min(distancedata_cluster_new$CoordinateZ)+
                 max(distancedata_cluster_new$CoordinateZ))
  center <- c(0,x_bar,y_bar,z_bar,0,0,0,0)
  M_new <- rbind(center, M_new)
    
  return(M_new)
  
}

plotcentermap <- function(center, M_new){
  Match <- matrix(rep(0,nrow(M_new)*2),ncol = 2)
  for (i in 1:nrow(M_new)){
    v <- solveequ(M_new, i)
    v_e <- create_expectV(center,M_new,i)
    distc2p <- center2point(center, M_new, i)
    anglediff <- as.numeric(findanglexy(v, v_e))
    #Match[i]<-anglediff
    Match[i,1] <-  anglediff #实际z轴与拟合法向夹角
    Match[i,2] <- distc2p  #实际rubisco与拟合羧酶体中心距离
  }
  return(Match)
}
