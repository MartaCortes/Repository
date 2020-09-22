setwd("~/Documents/Marta/MÁSTER DATA SCIENCE/Programming in R")
data=read.csv("iris_dataset.csv")


#My function
Mykmean <- function(x, y, Nclus, distance) {
  
#Random centroids assignament
  xcentroid <- runif(n=Nclus, min=min(x), max=max(x))
  ycentroid <- runif(n=Nclus, min=min(y), max=max(y))
  
#Storing of the centroids and coordinates info
  infocens <- data.frame(xcoor=x, ycoor=y, clusgroup = NA, error=NA)
  error <- data.frame(error=NULL)
  clusgroup <- data.frame(group = 1:Nclus, xcentroid = xcentroid, ycentroid=ycentroid)
  
  set <- FALSE
    loop=1
    while(set == FALSE) {
    
    # reasignamiento del grupo al centroide mas cercano
    for (i in 1:length(x)) {
      
#Distance calculation of each vector to the closest centroid (depending on the distance type)      
      if(distance=="Euclidean") {
         dist <- sqrt((x[i]-clusgroup$xcentroid)^2 + (y[i]-clusgroup$ycentroid)^2)
      }
         
      if(distance=="SEuclidean") {
        dist <- ((x[i]-clusgroup$xcentroid)^2 + (y[i]-clusgroup$ycentroid)^2)
      }
      
      if(distance=="Manhattan") {
        dist <- (abs(x[i]-clusgroup$xcentroid) + abs(y[i]-clusgroup$ycentroid))
      }
      
      if(distance=="Maximun") {
        dist <- max(abs(x[i]-clusgroup$xcentroid), abs(y[i]-clusgroup$ycentroid))
      }
      
#Error calculation of each iteration depending on the distance      
      infocens$clusgroup[i] <- which.min(dist)
      infocens$error[i] <- min(dist)
    }
      
    error=rbind(error,sum(infocens$error))    
    
#Centroids groups reassignament     
    xcentroid_past <- clusgroup$xcentroid
    ycentroid_past <- clusgroup$ycentroid
    
#Calculation of the new centers for each centroid 
    for(i in 1:Nclus) {
      clusgroup[i,2] <- mean(subset(infocens$xcoor, infocens$clusgroup == i))
      clusgroup[i,3] <- mean(subset(infocens$ycoor, infocens$clusgroup == i))
    }
    
#Adding of new data for the missing or "non-assigning" centroids
    if (length(unique(infocens$clusgroup))!=Nclus) {
      d <- Nclus-length(unique(infocens$clusgroup))
      for (i in 1:d) {
        addvars <- data.frame(group=which(is.na(clusgroup[,2])), xcentroid = runif(n=d, min=min(x), max=max(x)), 
                              ycentroid = runif(n=d, min=min(x), max=max(x)))
      }
      clusgroup <- rbind(clusgroup,addvars)
      clusgroup <- clusgroup[complete.cases(clusgroup), ]
    }
    
#Computing of the plot for each iteration
    plot=ggplot()+
      geom_point(infocens,mapping=aes(x=xcoor, y=ycoor, color=as.factor(clusgroup)))+
      geom_point(clusgroup, mapping=aes(x=xcentroid,y=ycentroid, color=as.factor(group)), shape=2, size=5)+
      labs(x="Axis x", y="Axis y")+
      theme(legend.title=element_blank())
    print(plot) 
    
#Printing of the number of iterations 
    loop=loop+1
    print(loop)
      
#Loop finalization if there is not changes in the centroid centers 
    if (identical(xcentroid_past, clusgroup$xcentroid) & identical(ycentroid_past, clusgroup$ycentroid)) set <- TRUE
    }
           
#Print of the plots and storing of the solution information (info about the centroids data, info about the centroids coordinates and error list)
   solution=list(Data=infocens, Clusters=clusgroup, Error=error)
   return(solution)
    
  
}

#Function implementation
a=Mykmean(data[,1],data[,3],5,"Euclidean")

