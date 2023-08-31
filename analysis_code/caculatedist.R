rm(list = ls())
PATH = getwd()
dir_fun = paste(PATH,"/analysis_code/",sep="")
dir = paste(dir_fun,"/data/",sep="")
# 计算距离的函数
calculate_distance <- function(df, i) {
  point1 <- df[i, c("x", "y", "z")]
  point2 <- df[i + 1, c("x", "y", "z")]
  distance <- sqrt(sum((point2 - point1)^2))
  return(distance)
}
data = read.table(paste(dir,"linker.txt",sep=""),head=FALSE)
dist = rep(0,2568)
colnames(data) = c("x", "y", "z")
for(i in 1:2568)
  if (i%%2==1){
    dist[i] = calculate_distance(data,i)
  }
    
data$dist = dist
write.table(data, paste(dir,"dist.txt",sep=""), sep = "\t", row.names = FALSE)
