
#SERIAL

# Introduction ------------------------------------------------------------

rm(list=ls())
setwd("~/201901_SDS_UC3M/3_SDC/Assignment1")

#install packages if required
required_packages <- c("dplyr",
                       "foreach",
                       "skimr",
                       "factoextra",
                       "cluster", 
                       "FactoMineR",
                       "psych",
                       "tidyr",
                       "caret",
                       "RColorBrewer", 
                       "parallel" ,
                       "doParallel") 

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

for (package in required_packages){
  library(package, character.only = T)
}

computers <-  read.csv2("computers.csv", row.names = 1,stringsAsFactors = F) 

computers$X <- NULL

#Change variables to binary
computers$multi <- ifelse(computers$multi == "no", 0,1)
computers$cd <- ifelse(computers$cd == "no", 0,1)
computers$premium <- ifelse(computers$premium == "no", 0,1)

#There are no NA

#Descriptive Analysis
head(computers)
skim(computers)

#Normalization needed as there are different measurement units
norm_data <- scale(computers) 


# K selection -------------------------------------------------------------


#Elbow Method for finding the optimal number of clusters

# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
kmeans_k <- function(k){
  set.seed(123456789)
  kmeans(norm_data, k, nstart=100,iter.max = 1000,)$tot.withinss
}

start.time <- Sys.time();
wss <- sapply(2:k.max, kmeans_k)
end.time <- Sys.time();
cat('Time to find the best value for k in serial programming is',(end.time - start.time),'minutes')

############### Elbow plot

ggplot(data.frame(wss),aes(x = 2:k.max, y = wss))+
  geom_path()+
  geom_point(shape=21, fill="white", color="black", size=3)+
  xlab("Number of clusters (K)")+
  ylab("Total within-clusters sum of squares")+
  scale_x_continuous(breaks = round(seq(2, k.max, by = 2),1))+
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",size = 2, linetype = "solid"))

#Confirm from The elbow graph - Check from 4 to 6
silhouette_score <- function(k){
  set.seed(123456789)
  km <- kmeans(norm_data, k, nstart=100,iter.max = 1000)
  ss <- silhouette(km$cluster, dist(norm_data))
  mean(ss[, 3])
}

start.time <- Sys.time();
avg_sil <- sapply(3:5, silhouette_score)
end.time <- Sys.time();
cat('Time to find the best value for k in serial programming is',(end.time - start.time),'minutes')

#Silhoutte plot

ggplot(data.frame(x = 3:5, avg_sil),aes(x = x, y = avg_sil))+
  geom_path()+
  geom_point(shape=21, fill="white", color="black", size=3)+
  xlab("Number of clusters (K)")+
  ylab("Average Silhouette Scores")+
  scale_x_continuous(breaks = round(seq(3, 5, by = 1),1))+
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",size = 2, linetype = "solid"))
 

#The maximum value is obtained when k = 4


# K-means algorithm -------------------------------------------------------

# Get K-means for K=4 and 100 initial random solutions and a maximum of 1000 iterations 
kmeans.data <- kmeans(norm_data,centers=4,iter.max=1000,nstart=100)

# The function returns a large list of objects. 

# The final partition is given in cluster
clusters <-  kmeans.data$cluster

# The sample mean vectors of the clusters are in centers
centers <- kmeans.data$centers

# Cluster plot ------------------------------------------------------------

#Plot the clusters using the 2 first PC

pca.data <-prcomp(norm_data)
get_eigenvalue(pca.data)

pca.data <- pca.data$x[,1:2]

data.cluster<- data.frame(pca.data)
data.cluster$cluster = as.factor(clusters)

ggplot(data.cluster, aes(x = PC1, y = PC2, fill = cluster, color=cluster))+
  geom_point(size = 3, shape = 21)+
  scale_fill_discrete(name = "Clusters")+
  scale_color_discrete(guide=FALSE)+
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",size = 2, linetype = "solid"))


# Mean Cluster ----------------------------------------------------------

computers$clusters <- clusters

mean_cluster<- function(k){
  comp_mean<- computers[clusters==k,]
  return(mean(comp_mean$price))
}

start.time <- Sys.time();
price <- c()
for(k in 1:4){
  price[k] <- mean_cluster(k)
}
end.time <- Sys.time();
cat('Time to find the average for each cluster is',(end.time - start.time),'seconds')


# HeadMap -----------------------------------------------------------------

df_kmeans= data.frame(cluster = c(1:4), centers)

norm_gather<- gather(df_kmeans, features, values, price:trend)

pal <-colorRampPalette(c("white","lightblue","navyblue"))

ggplot(data = norm_gather, aes(x = features, y = cluster, fill = values)) +
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +
  geom_tile() +
  xlab("Clusters")+
  ylab("Features")+
  scale_color_discrete(guide=FALSE)+
  scale_fill_gradientn(colours = pal(10),name = "Values")+
  theme(panel.background = element_rect(NULL))





#PARALLELIZATION

# Introduction ------------------------------------------------------------

rm(list=ls())
setwd("~/Downloads")

#install packages if required
required_packages <- c("dplyr",
                       "foreach",
                       "skimr",
                       "factoextra",
                       "cluster", 
                       "FactoMineR",
                       "psych",
                       "tidyr",
                       "caret",
                       "RColorBrewer", 
                       "parallel" ,
                       "doParallel") 

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

for (package in required_packages){
  library(package, character.only = T)
}

no_cores <- detectCores()

computers <-  read.csv2("computers.csv", row.names = 1,stringsAsFactors = F) 

computers$X <- NULL

#Change variables to binary
computers$multi <- ifelse(computers$multi == "no", 0,1)
computers$cd <- ifelse(computers$cd == "no", 0,1)
computers$premium <- ifelse(computers$premium == "no", 0,1)

#There are no NA

#Descriptive Analysis
head(computers)
skim(computers)

#Normalization needed as there are different measurement units
norm_data <- scale(computers) 


# K selection -------------------------------------------------------------


#Elbow Method for finding the optimal number of clusters

# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
kmeans_k <- function(k){
  set.seed(123456789)
  kmeans(norm_data, k, nstart=100,iter.max = 1000,)$tot.withinss
}

start.time <- Sys.time();
clust <- makeCluster(no_cores, type='PSOCK')
clusterExport(clust, "norm_data")
wss=parSapply(cl=clust,2:15, FUN=kmeans_k)
stopCluster(clust)
end.time <- Sys.time();
cat('Time to find the best value for k in cluster programming using psock and parSapply is',difftime(end.time, start.time, units = "secs"),'seconds')

start.time <- Sys.time();
registerDoParallel(makeCluster(no_cores,type = "PSOCK"))
wss=foreach(i=2:15,.combine='cbind')%dopar%{
  kmeans_k(i)
}
stopImplicitCluster()
end.time <- Sys.time();
cat('Time to find the best value for k in cluster programming using psock and foreach is',difftime(end.time, start.time, units = "secs"),'seconds')

start.time <- Sys.time();
clust <- makeCluster(no_cores, type='FORK')
wss=parSapply(cl=clust,2:15, FUN=kmeans_k)
stopCluster(clust)
end.time <- Sys.time();
cat('Time to find the best value for k in cluster programming using fork and parSapply is',difftime(end.time, start.time, units = "secs"),'seconds')

start.time <- Sys.time();
registerDoParallel(makeCluster(no_cores,type = "FORK"))
wss=foreach(i=2:15,.combine='cbind')%dopar%{
  kmeans_k(i)
}
stopImplicitCluster()
end.time <- Sys.time();
cat('Time to find the best value for k in cluster programming using psock and foreach is',difftime(end.time, start.time, units = "secs"),'seconds')


############### Elbow plot

ggplot(data.frame(wss),aes(x = 2:k.max, y = wss))+
  geom_path()+
  geom_point(shape=21, fill="white", color="black", size=3)+
  xlab("Number of clusters (K)")+
  ylab("Total within-clusters sum of squares")+
  scale_x_continuous(breaks = round(seq(2, k.max, by = 2),1))+
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",size = 2, linetype = "solid"))

#The maximum value is obtained when k = 4


# K-means algorithm -------------------------------------------------------

# Get K-means for K=4 and 100 initial random solutions and a maximum of 1000 iterations 
kmeans.data <- kmeans(norm_data,centers=4,iter.max=1000,nstart=100)

# The function returns a large list of objects. 

# The final partition is given in cluster
clusters <-  kmeans.data$cluster

# The sample mean vectors of the clusters are in centers
centers <- kmeans.data$centers

# Cluster plot ------------------------------------------------------------

#Plot the clusters using the 2 first PC

pca.data <-prcomp(norm_data)
get_eigenvalue(pca.data)

pca.data <- pca.data$x[,1:2]

data.cluster<- data.frame(pca.data)
data.cluster$cluster = as.factor(clusters)

ggplot(data.cluster, aes(x = PC1, y = PC2, fill = cluster, color=cluster))+
  geom_point(size = 3, shape = 21)+
  scale_fill_discrete(name = "Clusters")+
  scale_color_discrete(guide=FALSE)+
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",size = 2, linetype = "solid"))


# Mean Cluster ----------------------------------------------------------

computers$clusters <- clusters

mean_cluster<- function(k){
  comp_mean<- computers[computers$clusters==k,]
  return(mean(comp_mean$price))
}

start.time <- Sys.time();
registerDoParallel(makeCluster(no_cores,type = "PSOCK"))
price <- foreach(i=1:4,.combine='cbind')%dopar%{
  mean_cluster(i)
}
stopImplicitCluster()
end.time <- Sys.time();
cat('Time to find the average for each cluster using PSCOCK is',difftime(end.time, start.time, units = "secs"),'seconds')

start.time <- Sys.time();
registerDoParallel(makeCluster(no_cores,type = "FORK"))
price <- foreach(i=1:4,.combine='cbind')%dopar%{
  mean_cluster(i)
}
stopImplicitCluster()
end.time <- Sys.time();
cat('Time to find the average for each cluster using FORK is',difftime(end.time, start.time, units = "secs"),'seconds')


# HeadMap -----------------------------------------------------------------

df_kmeans= data.frame(cluster = c(1:4), centers)

norm_gather<- gather(df_kmeans, features, values, price:trend)

pal <-colorRampPalette(c("white","lightblue","navyblue"))

ggplot(data = norm_gather, aes(x = features, y = cluster, fill = values)) +
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +
  geom_tile() +
  xlab("Clusters")+
  ylab("Features")+
  scale_color_discrete(guide=FALSE)+
  scale_fill_gradientn(colours = pal(10),name = "Values")+
  theme(panel.background = element_rect(NULL))





