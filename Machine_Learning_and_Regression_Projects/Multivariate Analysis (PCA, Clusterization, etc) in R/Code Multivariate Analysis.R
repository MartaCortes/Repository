
rm(list=ls())

setwd("~/Documents/Marta/MÁSTER DATA SCIENCE/Multivariate Analysis/Assignment 1 Multivariate/")
library(readxl)
library(xlsx)
library(dplyr)
library("ggplot2")
library("GGally")
library("psych")
library("reshape2")
library("ISLR")
library("factoextra")
library("plotrix")
library(MASS)
library("corrplot")
library(skimr)

datat <- read.csv2("Data set ACS.csv")
data_region <- read_excel("States per region.xlsx")


###### Treatment for the data set
#Impute 0 to NA case
datat[is.na(datat)]=0
#datat=mutate(datat,Women=Women/TotalPop*100, VotingCitizen=VotingCitizen/TotalPop*100)
names_subset=c("County","State","Women","Hispanic", "White", "Black", "VotingCitizen","IncomePerCap","Poverty","Professional","Service","Office","Construction","Drive","Transit","Walk","MeanCommute","Unemployment","PublicWork","PrivateWork",'SelfEmployed')
data=select_(datat, .dots=names_subset)

variables=ncol(data)
observations=nrow(data)

summary(data)
skim(data)

#individual analysis of the variables

#Plots of densities of all variables in order to see the distributions
summarydata_melt <- melt(data[,-c(1,2)])

ggplot(summarydata_melt, aes (value, fill=variable)) +
  geom_density() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")


#Plots QQplot
library(fitdistrplus)
ggplot(summarydata_melt, aes (sample=value, color=variable)) +
  geom_qq() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")

#Plots Boxplots
ggplot(data = summarydata_melt, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() +
  facet_wrap( ~ variable, scales="free") +
  theme(legend.position="none")


# Parallel coordinates plots and andrews
datar=left_join(data,data_region)
datar$Region=as.factor(datar$Region)
colors <- c("black","darkblue","green","cyan","purple")[as.factor(datar$Region)]
parcoord(datar[,-c(1,2,22)],col=colors,var.label=TRUE)


ggparcoord(datar,
           columns = 3:21, groupColumn = 22,
           showPoints = TRUE, 
           alphaLines = 0.3,
           scale="center"
) + 
  theme(
    plot.title = element_text(size=10)
  )+theme(axis.text.x = element_text(angle = 90)) 

library("pracma")
andrewsplot(as.matrix(datar[,-c(1,2,22)]),datar$Region,style="cart",npts=3220)

#library(andrews)
#dataa=datar[,-c(1,2)]
#dataa$Region=as.character(dataa$Region)
#andrews(as.matrix(dataa),type=1,clr=20)

#Descriptive
library(skimr)
skim(data)
summary(data)

#Transformations

library(rcompanion)
library(EnvStats)
transformTukey(100-data$White, plotit = FALSE)
white2=((100-data$White)^(0.175))
qqPlot(white2, add.line = TRUE)

transformTukey((10+data$Black), plotit = FALSE)
black2=-1*(10+data$Black)^(-3.125)
qqPlot(black2, add.line = TRUE)

transformTukey((1+data$Transit), plotit = FALSE)
transit2=-1*(1+data$Transit)^(-1.65)
qqPlot(transit2, add.line = TRUE)

data1=mutate(data, Women=log(Women),Hispanic=log(Hispanic+1), White=(100-White)^0.175, Black=(Black)^0.225,
             VotingCitizen=log(VotingCitizen), Poverty=log(Poverty),
             Construction=log(Construction+1),Drive=(Drive)^3,Transit=(-1*(1+Transit)^(-1.65)),Walk=log(Walk+1),Unemployment=log(Unemployment+1),
             PublicWork=log(PublicWork),PrivateWork=(PrivateWork)^3,SelfEmployed= log(SelfEmployed+1))


#Plots after transformations
summarydata_melt <- melt(data1[,-c(1,2)])

ggplot(summarydata_melt, aes (value, fill=variable)) +
  geom_density() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")

ggplot(summarydata_melt, aes (sample=value, color=variable)) +
  geom_qq() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")

#Plot Boxplots
ggplot(summarydata_melt, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() +
  facet_wrap( ~ variable, scales="free") +
  theme(legend.position="none")



#detecting outliers
#Squared Mahalanobis distances with respect to the sample mean vector for counties
datan=mutate(data1, Name=paste0(County,"-",State))
datan_names=datan$Name
datac=datan[,-c(1,2)]
datac=as.matrix(datac)
datac=apply(datac,2,as.numeric)
datac=datac[,-20]
rownames(datac)=datan_names

require("robustbase")

##################################################
#Squared Mahalanobis distances with respect to the sample mean vector for the counties

#Minimum Covariance Determinant (MCD) estimators

mcd.datac <- covMcd(datac,alpha=0.90)

rob.m.datac <- mcd.datac$center
rob.m.datac

rob.S.datac <- mcd.datac$cov
rob.S.datac

#Robust Mahalanobis distances

rob.mah.datac <- mahalanobis(datac,rob.m.datac,rob.S.datac)
rob.mah.datac
sort.rob.mah.datac <- sort(rob.mah.datac,index.return=TRUE)$x
sort.rob.mah.datac

plot(sort.rob.mah.datac,pch=19,col="blue",xlab="",ylab="",main="Robust Mahalanobis distances for counties")

# FDA
n.datac <- nrow(datac)
n.datac
p.datac <- ncol(datac)
p.datac

p.values.rob.datac <- 1 - pchisq(rob.mah.datac,p.datac)
p.values.rob.datac

# Sort them in increasing order

sort.p.values.datac <- sort(p.values.rob.datac,index.return=TRUE)
#Determination of outliers

which(sort.p.values.datac$x < ((1:n.datac)/n.datac*0.01))

mah.rob.countiesout=as.data.frame(which(sort.p.values.datac$x < ((1:n.datac)/n.datac*0.01)))

mah.rob.countiesout$Name=rownames(mah.rob.countiesout)

library(tidyr)
#mah.rob.countiesout=separate(mah.rob.countiesout,into="County",Name, sep="-", remove=FALSE,convert=TRUE)
mah.rob.countiesout=left_join(mah.rob.countiesout, dplyr::select(datan, County, State,Name))


#Counties by state
data_states=datat%>%group_by(State)%>%summarise(Counties=n())
countiesout_states=mah.rob.countiesout%>%group_by(State)%>%summarise(Counties_out=n())
data_states=left_join(data_states,countiesout_states)
data_states=mutate(data_states,Percentage_out=Counties_out/Counties)

data_states2=left_join(data_states, data_region)
data_region2=data_states2%>%group_by(Region)%>%summarise(Counties=sum(Counties, na.rm=TRUE),Counties_out=sum(Counties_out, na.rm=TRUE), Percentage_out=sum(Percentage_out, na.rm=TRUE))
#eliminate the outliers before estimating the covariance and correlation matrices
#data1=data1[!data1$State%in%c("Hawaii","Alaska","District of Columbia","Puerto Rico"),]
data1=data1[!data1$State%in%c("Puerto Rico"),]

#### data1 final data ya transformada y quitada los outliers escogidos

#Standarizations
scaled.data1 <- scale(data1[,3:21])
scaled.data1 <- as.data.frame(scaled.data1)
scaled.data2 <- cbind(data1[,1:2],scaled.data1)

#Robust estandarization
#library(quantable)
#scaled.data1=robustscale(data1[,-c(1,2)], dim = 2, center = TRUE, scale = TRUE,
#            preserveScale = FALSE)


#Check that we get mean of 0 and sd of 1
#colMeans(scaled.data1$data) 
#Faster version of apply(scaled.dat, 2, mean)
#apply(scaled.data1$data, 2, sd)


#Check that we get mean of 0 and sd of 1
colMeans(scaled.data1) 
#Faster version of apply(scaled.dat, 2, mean)
apply(scaled.data1, 2, sd)


#Correlation Plots of all variables (data data without outliers) ------------------------------>>>>>>>>>>>>>>>>>>>>
pairs.panels(scaled.data1[,1:19], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             cex=5,
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)



####################
#PCA ##############
####################

#################################################################################################################
# Sample sizes and dimensions of the scale data set

X.PCA <- scaled.data1

n.X.PCA <- nrow(X.PCA)
n.X.PCA
p.X.PCA <- ncol(X.PCA)
p.X.PCA

##################################################################################################################
# PCA for the USA Census data set
##################################################################################################################

##################################################################################################################
# There are several functions in R to perform PCA. Next, we use the function prcomp

PCS.X.PCA <- prcomp(X.PCA)

# Have a look at the outputs

names(X.PCA)
summary(PCS.X.PCA)
##################################################################################################################
# The PC scores is are in x. In this case, the scores are a matrix of size 19x19

dim(PCS.X.PCA$x)
head(PCS.X.PCA$x)

##################################################################################################################
# Make a plot of the first PC scores

plot(PCS.X.PCA$x[,1],pch=19,col="deepskyblue2",main="First PCA")
abline(h=0)

# Make a plot of the second PC scores
plot(PCS.X.PCA$x[,2],pch=19,col="deepskyblue2",main="Second PCs")
abline(h=0)

# Make a plot of the three PC scores
plot(PCS.X.PCA$x[,3],pch=19,col="deepskyblue2",main="Third PCs")
abline(h=0)

# Make a plot of the four PC scores
plot(PCS.X.PCA$x[,4],pch=19,col="deepskyblue2",main="Fourth PCs")
abline(h=0)

# Make a plot of the five PC scores
plot(PCS.X.PCA$x[,5],pch=19,col="deepskyblue2",main="Fifth PCs")
abline(h=0)

##################################################################################################################
# Make a plot of the first two PC scores
plot(PCS.X.PCA$x[,1:2],pch=19,col="deepskyblue2",main="First two PCs")
abline(h=0,v=0)

# Make a plot of the second and third PC scores
plot(PCS.X.PCA$x[,2:3],pch=19,col="deepskyblue2",main="Second and third two PCs")
abline(h=0,v=0)

# Make a plot of the third and fourth PC scores
plot(PCS.X.PCA$x[,3:4],pch=19,col="deepskyblue2",main="Third and Fourth two PCs")
abline(h=0,v=0)

# Make a plot of the fourth and fifth PC scores
plot(PCS.X.PCA$x[,4:5],pch=19,col="deepskyblue2",main="Fourth and fifth two PCs")
abline(h=0,v=0)

# Make a plot of the first and third PC scores
plot(PCS.X.PCA$x[,1],PCS.X.PCA$x[,3],pch=19,col="deepskyblue2",main="First and third two PCs")
abline(h=0,v=0)

# Make a plot of the second and fourth PC scores
plot(PCS.X.PCA$x[,2],PCS.X.PCA$x[,4],pch=19,col="deepskyblue2",main="Second and fourth two PCs")
abline(h=0,v=0)

# Make a plot of the first and fourth PC scores
plot(PCS.X.PCA$x[,1],PCS.X.PCA$x[,4],pch=19,col="deepskyblue2",main="First and fourth two PCs")
abline(h=0,v=0)

# Make a plot of the first and fifth PC scores
plot(PCS.X.PCA$x[,1],PCS.X.PCA$x[,5],pch=19,col="deepskyblue2",main="First and fifth two PCs")
abline(h=0,v=0)

# Make a plot of the second and fifth PC scores
plot(PCS.X.PCA$x[,2],PCS.X.PCA$x[,5],pch=19,col="deepskyblue2",main="Second and fifth two PCs")
abline(h=0,v=0)



install.packages("devtools")
install.packages("usethis")
library(devtools)
install_github("vqv/ggbiplot")

install.packages("ggbiplot")
library(ggbiplot)
library(grid)

#Plot of first two PC scores with the initial variables
biplot(PCS.X.PCA,col=c("deepskyblue2","firebrick2"),cex=c(0.5,0.8))

##################################################################################################################
# The eigenvalues of the sample covariance matrix of X, i.e., the variances of the PCs are the square of sdev
# As in this example n<p, only 5 eigenvalues can be different from 0. These are those that appear here.

PCS.X.PCA$sdev^2

##################################################################################################################
# Have a look at these eigenvalues

# Screeplot with the 19 eigenvalues

fviz_eig(PCS.X.PCA,ncp=19,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")

# Screeplot with the first 5 eigenvalues

fviz_eig(PCS.X.PCA,ncp=5,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")

##################################################################################################################
# How many PCs are important?

# Have a look at the proportion of explained variance and the cumulative proportion of explained variance

get_eigenvalue(PCS.X.PCA)

# We need 6 PCs to have the 75% of the total variability of the data set
# Have a look at the sample mean of the eigenvalues

eval.PCS.X.PCA <- PCS.X.PCA$sdev^2
mean(eval.PCS.X.PCA)

# The number of eigenvalues larger than this sample mean is 

sum(eval.PCS.X.PCA>mean(eval.PCS.X.PCA))

# Then, the dimension of the data set is reduced from 19 to either 5 that represents the number of variables 
# in the original data set keeping around the 70% of the information inside

##################################################################################################################
# The loading matrix, i.e., the eigenvectors of the sample covariance matrix of X are given in rotation

dim(PCS.X.PCA$rotation)
PCS.X.PCA$rotation[,1:5]

# Note that only 19 eigenvectors (corresponding to 19 PCs) appear

##################################################################################################################
# Interpretation of the PCs - Hacer con los PCAs que más diferentes pesos tengan

# Each value represent the weight of the associated variable
# Have a look at the important variables in the first PC

plot(PCS.X.PCA$rotation[,1],pch=20,col="deepskyblue2",main="Weights for the first PC")
abline(h=0)

# Have a look at the important variables in the second PC

plot(PCS.X.PCA$rotation[,2],pch=20,col="deepskyblue2",main="Weights for the second PC")
abline(h=0)

# Have a look at the important variables in the first two PCs

plot(PCS.X.PCA$rotation[,1:2],pch=19,col="deepskyblue2",main="Relevant Variables for the first two PCs")
abline(h=0,v=0)
text(PCS.X.PCA$rotation[,1:2],labels=colnames(X.PCA),pos = 1,col="firebrick2",cex=0.5)

draw.circle(0,0,0.29,border="green2",lwd=4)


## Plot the scores of the five PCs for a qualitative variable (Region)
# The PCs will show that the variables have different behavior in terms of the four groups

scaled.data3 <- scaled.data2
scaled.data3[,2] <- as.character(scaled.data3[,2])

a <- c("Alabama", "Arkansas", "Delaware", "District of Columbia", "Florida", "Georgia", "Kentucky", "Louisiana", "Maryland", 
       "Mississippi", "North Carolina", "Oklahoma", "South Carolina", "Tennessee", "Texas", "Virginia", "West Virginia")

b <- c("Illinois", "Indiana", "Iowa", "Kansas", "Michigan", "Minnesota", "Missouri", "Nebraska", "North Dakota", "Ohio","South Dakota", "Wisconsin")

c <- c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "New Jersey", "New York", "Pennsylvania", "Rhode Island", "Vermont")

d <- c("Alaska", "Arizona", "California", "Colorado", "Hawaii", "Idaho", "Montana", "Nevada", "New Mexico", "Oregon", "Utah", "Washington","Wyoming")


for (i in 1:length(scaled.data3[,2])) {
  if (((scaled.data3[i,2]) %in% a) == TRUE) {
    scaled.data3[i,2] <- "South"
  }
  if (((scaled.data3[i,2]) %in% b) == TRUE) {
    scaled.data3[i,2] <- "Midwest"
  }
  if (((scaled.data3[i,2]) %in% c) == TRUE) {
    scaled.data3[i,2] <- "Northeast"
  }
  if (((scaled.data3[i,2]) %in% d) == TRUE) {
    scaled.data3[i,2] <- "West"
  }
}

scaled.data3[,2] <- as.factor(scaled.data3[,2])


Y <- scaled.data3[,2]
colors.X <- c("royalblue1","firebrick1","forestgreen","darkorchid")[Y]
parcoord(X.PCA,col=colors.X,var.label=TRUE)

pairs(X.PCA,pch=20,col=colors.X)
pairs(PCS.X.PCA$x[,1:5],col=colors.X,pch=19,main="The first five PCs", oma=c(3,3,3,15))
par(xpd = TRUE)
legend("right", fill = unique(c("royalblue1","firebrick1","forestgreen","darkorchid")), legend = levels(Y),cex=0.5)


##################################################################################################################
# Plot the correlations between the original data set and the scores
# This is useful to understand also which are the important variables in the data set in terms of variability

corrplot(cor(X.PCA,PCS.X.PCA$x),is.corr=T)

# Reduce to only five PCs

corrplot(cor(X.PCA,PCS.X.PCA$x[,1:5]),is.corr=T)

##################################################################################################################
#CLUSTER ANALYSIS
##################################################################################################################

# Load the data set into memory
scaled.data1

# Define the data matrix

X.scaled.data1 <- scaled.data1

# Sample size and dimension of the Census US data set

n.X.scaled.data1 <- nrow(X.scaled.data1)
n.X.scaled.data1
p.X.scaled.data1 <- ncol(X.scaled.data1)
p.X.scaled.data1

##################################################################################################################
# Remember: PCA for the Census US data set

PCS.X.scaled.data1 <- prcomp(X.scaled.data1)

##################################################################################################################
# Make a plot of the first two PCs

par(mfrow=c(1,1))
plot(PCS.X.scaled.data1$x[,1:2],pch=19,col="deepskyblue2")

# This plot suggests that there migth be some heterogeneity, i.e., at least two clusters

# Remember that around 5 PCs where needed to explain at least the 70% of the total variability of the data set
# Thus, the number of clusters might be more than 2

##################################################################################################################
# Try K=2 and 100 initial random solutions

kmeans.X.scaled.data1 <- kmeans(X.scaled.data1,centers=2,iter.max=1000,nstart=100)

# The function returns a large list of objects. 
summary(kmeans.X.scaled.data1)

# The final partition is given in cluster
kmeans.X.scaled.data1$cluster

# The number of observations in each cluster
kmeans.X.scaled.data1$size

# The sample mean vectors of the clusters are in centers
kmeans.X.scaled.data1$centers

# Total within-cluster sum of squares
kmeans.X.scaled.data1$tot.withinss 

# Make a plot of the first two PCs split in these two clusters
colors.kmeans.X.scaled.data1 <- c("deepskyblue2","firebrick2")[kmeans.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.kmeans.X.scaled.data1)

# It does not appear to be a bad solution but the limits of the clusters are a bit mixed

##################################################################################################################
# Compare several values of K and select the most appropriate one with #### WSS ####
fviz_nbclust(X.scaled.data1,kmeans,method="wss", k.max=10) +
  geom_vline(xintercept = 4, linetype = 2)

# The plot have the bend (knee) starting in the number 4, so it suggests the presence of 4 cluster

##################################################################################################################
# Compare several values of K and select the most appropriate one with #### Silhouette ####
fviz_nbclust(X.scaled.data1,kmeans,method="elbow",k.max=10)

# The plot suggests the presence of 2 clusters

##################################################################################################################
# Compare several values of K and select the most appropriate one with the #### gap statistic ####
fviz_nbclust(X.scaled.data1,kmeans,method="gap",k.max=10,nboot=100) +
  geom_vline(xintercept = 4, linetype = 2)

# The plot suggests the presence of 10 clusters or 4 clusters


##################################################################################################################
# Let see the solution for K=2
kmeans2.X.scaled.data1 <- kmeans(X.scaled.data1,centers=2,iter.max=1000,nstart=100)

# The cluster solution
kmeans2.X.scaled.data1$cluster

# Make a plot of the first two PCs split in these two clusters
colors.kmeans2.X <- c("firebrick2", "deepskyblue2")[kmeans2.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.kmeans2.X)

# Another way to do the plot
fviz_cluster(kmeans2.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

##################################################################################################################
# Let see the solution for K=4
kmeans4.X.scaled.data1 <- kmeans(X.scaled.data1,centers=4,iter.max=1000,nstart=100)

# The cluster solution
kmeans4.X.scaled.data1$cluster

# Make a plot of the first two PCs split in these four clusters
colors.kmeans4.X <- c("deepskyblue2","forestgreen","firebrick2","orange")[kmeans4.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.kmeans4.X)

# Another way to do the plot
fviz_cluster(kmeans4.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

##################################################################################################################
# Let see the solution for K=3
kmeans3.X.scaled.data1 <- kmeans(X.scaled.data1,centers=3,iter.max=1000,nstart=100)

# The cluster solution
kmeans3.X.scaled.data1$cluster

# Make a plot of the first two PCs split in these tree clusters
colors.kmeans3.X <- c("firebrick2","forestgreen","deepskyblue2")[kmeans3.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.kmeans3.X)

# Another way to do the plot
fviz_cluster(kmeans3.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])


##################################################################################################################
# Silhouette plot for the solution
# Compute the silhouette

library("cluster")

#K=2
sil.kmeans2.X.scaled.data1 <- silhouette(kmeans2.X.scaled.data1$cluster,dist(X.scaled.data1,"euclidean"))
start.time1 <- Sys.time()
fviz_silhouette(sil.kmeans2.X.scaled.data1)
end.time1 <- Sys.time()
computation.time1 <-  end.time1 - start.time1

#Average silhouette width --> 1=0.22 and 2=0.12 (mean=0.19)

#K=4
sil.kmeans4.X.scaled.data1 <- silhouette(kmeans4.X.scaled.data1$cluster,dist(X.scaled.data1,"euclidean"))
fviz_silhouette(sil.kmeans4.X.scaled.data1)
#Average silhouette width --> 1= 0.22, 2=0.06, 3=0.15, 4=0.09 (mean=0.14)

#K=3
sil.kmeans3.X.scaled.data1 <- silhouette(kmeans3.X.scaled.data1$cluster,dist(X.scaled.data1,"euclidean"))
fviz_silhouette(sil.kmeans3.X.scaled.data1)
#Average silhouette width --> 1=0.12, 2=017, 3=0,17 (mean=0.16)

#Looks like the K=2 is the best (highest Average silhouette width)


##################################################################################################################
# Perform K-Means clustering for the principal components of the Census US data set
##################################################################################################################

# We take 5 PCs as in Topic 3 and K=2
kmeans.PCS.X.scaled.data1 <- kmeans(PCS.X.scaled.data1$x[,1:5],centers=2,iter.max=1000,nstart=100)

# Make a plot of the first two PCs split in these two clusters
colors.kmeans.PCS.X.scaled.data1 <- c("firebrick2", "deepskyblue2")[kmeans.PCS.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.kmeans.PCS.X.scaled.data1)

# Another way to do the plot
fviz_cluster(kmeans.PCS.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

# Compute the silhouette
sil.PCS.X.scaled.data1 <- silhouette(kmeans.PCS.X.scaled.data1$cluster,dist(X.scaled.data1,"euclidean"))
fviz_silhouette(sil.PCS.X.scaled.data1)
#Average silhouette width --> red=0.12 and blue=0.22 (mean=0.19)

# Thus, the solution given with the PCs is very close to the one with the original data
# Note that the original data set uses 19 variables while the PCs only uses 5 variables!!!

##################################################################################################################
##################################################################################################################
# Perform K-Medoids clustering for the Census US data set
##################################################################################################################
##################################################################################################################

##################################################################################################################
# Determine the optimal k by optimum average silhouette width for medoids

library(fpc)
pamk.best <- pamk(X.scaled.data1)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n") +
  title("number of clusters estimated by optimum average silhouette width", line = -2)
plot(pam(X.scaled.data1, pamk.best$nc))

#The plot suggest k=3

##################################################################################################################
# Run PAM with Euclidean distance for k=2 and k=3
##################################################################################################################
#K=2
pam2.E.X.scaled.data1 <- pam(X.scaled.data1,k=2,metric="euclidean",stand=FALSE)

# Make a plot of the first two PCs split in these two clusters
colors.pam2.E.X.scaled.data1 <- c("deepskyblue2","firebrick2")[pam2.E.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.pam2.E.X.scaled.data1)

# Another way to do the plot
fviz_cluster(pam2.E.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

# Have a look at the silhouette
sil.pam2.E.X.scaled.data1<- silhouette(pam2.E.X.scaled.data1$cluster,dist(X.scaled.data1,method="euclidean"))
fviz_silhouette(sil.pam2.E.X.scaled.data1)

# The solution is sligthly worse than the one with K-means (Average Silhouette width = 0.12)

#K=3
pam3.E.X.scaled.data1 <- pam(X.scaled.data1,k=3,metric="euclidean",stand=FALSE)

# Make a plot of the first two PCs split in these two clusters
colors.pam3.E.X.scaled.data1 <- c("deepskyblue2","firebrick2", "forestgreen")[pam3.E.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.pam3.E.X.scaled.data1)

# Another way to do the plot
fviz_cluster(pam3.E.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

# Have a look at the silhouette
sil.pam3.E.X.scaled.data1<- silhouette(pam3.E.X.scaled.data1$cluster,dist(X.scaled.data1,method="euclidean"))
fviz_silhouette(sil.pam3.E.X.scaled.data1)

# The solution is sligthly worse than the one with K-means (Average Silhouette width = 0.14) but better than the one with k=2

##################################################################################################################
# Run PAM with Manhattan distance for k=2 and k=3
##################################################################################################################

#K=2
pam2.M.X.scaled.data1 <- pam(X.scaled.data1,k=2,metric="manhattan",stand=FALSE)

# Make a plot of the first two PCs split in these two clusters
colors.pam2.M.X.scaled.data1 <- c("deepskyblue2","firebrick2")[pam2.M.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.pam2.M.X.scaled.data1)

# Another way to do the plot
fviz_cluster(pam2.M.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

# Have a look at the silhouette
sil.pam2.M.X.scaled.data1<- silhouette(pam2.M.X.scaled.data1$cluster,dist(X.scaled.data1,method="manhattan"))
fviz_silhouette(sil.pam2.M.X.scaled.data1)

# The solution is sligthly better than the one with Euclidean (Average Silhouette width = 0.13) but worse than with Kmeans

#K=3
pam3.M.X.scaled.data1 <- pam(X.scaled.data1,k=3,metric="manhattan",stand=FALSE)

# Make a plot of the first two PCs split in these two clusters
colors.pam3.M.X.scaled.data1 <- c("deepskyblue2","firebrick2", "forestgreen")[pam3.M.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.pam3.M.X.scaled.data1)

# Another way to do the plot
fviz_cluster(pam3.M.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

# Have a look at the silhouette
sil.pam3.M.X.scaled.data1<- silhouette(pam3.M.X.scaled.data1$cluster,dist(X.scaled.data1,method="manhattan"))
fviz_silhouette(sil.pam3.M.X.scaled.data1)

# The solution is the same than the one with Euclidean (ASW=0.14) but worse than with Kmeans


##################################################################################################################
# Run CLARA clustering for the Census US data set
##################################################################################################################

pamk.best <- pamk(X.scaled.data1, usepam=FALSE)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n") +
  title("number of clusters estimated by optimum average silhouette width", line = -2)
plot(pam(X.scaled.data1, pamk.best$nc))


#number of clusters estimated by optimum average silhouette width: 7

##################################################################################################################
# Fix K=2 and k=7 clusters
# K=2
clara2.X.scaled.data1 <- clara(X.scaled.data1,k=2,metric="manhattan",stand=FALSE,samples=100,sampsize=32)

# Make a plot of the first two PCs split in these two clusters

colors.clara2.X.scaled.data1 <- c("deepskyblue2","firebrick2")[clara2.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.clara2.X.scaled.data1)

# Another way to do the plot
fviz_cluster(clara2.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

# Have a look at the silhouette
sil.clara2.X.scaled.data1 <- silhouette(clara2.X.scaled.data1$cluster,dist(X.scaled.data1,method="manhattan"))
start.time2 <- Sys.time()
fviz_silhouette(sil.clara2.X.scaled.data1)
end.time2 <- Sys.time()
computation.time2 <-  end.time2 - start.time2
difference.time <- computation.time2 - computation.time1

# The solution is the same than the one with Kmeans (ASW=0.19) but the computational effort is higher (0.08500814 secs)

#K=7
clara7.X.scaled.data1 <- clara(X.scaled.data1,k=7,metric="manhattan",stand=FALSE,samples=100,sampsize=32)

# Make a plot of the first two PCs split in these two clusters

colors.clara7.X.scaled.data1 <- c("deepskyblue2","firebrick2", "forestgreen", "orange", "black", "pink", "darkmagenta")[clara7.X.scaled.data1$cluster]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.clara7.X.scaled.data1)

# Another way to do the plot
fviz_cluster(clara7.X.scaled.data1,data=PCS.X.scaled.data1$x[,1:2])

# Have a look at the silhouette
sil.clara7.X.scaled.data1 <- silhouette(clara7.X.scaled.data1$cluster,dist(X.scaled.data1,method="manhattan"))
fviz_silhouette(sil.clara7.X.scaled.data1)

# The solution is much worse than the one with k=2 (ASW=0.12)

##################################################################################################################
# Run K-Medoids clustering for mixed type data sets with the Census US data set (mixed data)
##################################################################################################################

##################################################################################################################
# Compute the Gower distance for the observations in the data set

dist.Gower.scaled.data2 <- daisy(scaled.data2,metric="gower")
summary(dist.Gower.scaled.data2)

# In the output, "I" means quantitative variable while N means qualitative variable

# Find the closest counties with the Gower distance
mat.dist.Gower.scaled.data2 <- as.matrix(dist.Gower.scaled.data2)
scaled.data2[which(mat.dist.Gower.scaled.data2==min(mat.dist.Gower.scaled.data2[mat.dist.Gower.scaled.data2!=min(mat.dist.Gower.scaled.data2)]),arr.ind = TRUE)[1,],]

# Find the most distant counties with the Gower distance
scaled.data2[which(mat.dist.Gower.scaled.data2==max(mat.dist.Gower.scaled.data2[mat.dist.Gower.scaled.data2!=max(mat.dist.Gower.scaled.data2)]),arr.ind = TRUE)[1,],]

##################################################################################################################
# Consider K=2:10, run PAM and select K using the silhouette

K.scaled.data2 <- matrix(NA,nrow=1,ncol=9)
for (i in 1:9){
  kmedoids.scaled.data2 <- pam(mat.dist.Gower.scaled.data2,k=i+1,diss=TRUE)
  K.scaled.data2[i] <- kmedoids.scaled.data2$silinfo$avg.width
}
plot(2:10,K.scaled.data2,pch=20,col="deepskyblue2")
which.max(K.scaled.data2)+1

#it suggests to have K=6

##################################################################################################################
# Run the algorithm for K=6 and get some information from the results

kmedoids.scaled.data2 <- pam(mat.dist.Gower.scaled.data2,k=6,diss=TRUE)

# Medoids
scaled.data2[kmedoids.scaled.data2$medoids,]

# Have a look at the silhouette
sil.kmedoids.scaled.data2 <- silhouette(kmedoids.scaled.data2$cluster,mat.dist.Gower.scaled.data2)
fviz_silhouette(sil.kmedoids.scaled.data2)
summary(sil.kmedoids.scaled.data2)

#Worst result than Kmeans (ASW=0.08)

##################################################################################################################
##################################################################################################################
# Hierarchical clustering analysis for the Census US data set
##################################################################################################################
##################################################################################################################

##################################################################################################################
# Agglomerative hierarchical clustering analysis for the Census US data set
##################################################################################################################

##################################################################################################################
# Compute the Euclidean distance matrix between the observations in the data matrix

dist.X.scaled.data1 <- daisy(X.scaled.data1,metric="euclidean",stand=FALSE)

##################################################################################################################
# Single linkage
single.X.scaled.data1 <- hclust(dist.X.scaled.data1,method="single")

# Plot dendogram of the solution and take k=2 as with Kmeans
plot(single.X.scaled.data1,main="Single linkage",cex=0.8)
rect.hclust(single.X.scaled.data1,k=2,border="deepskyblue2")


# See the assignment
cutree(single.X.scaled.data1,2)
table(cutree(single.X.scaled.data1,2))

# Make a plot of the first two PCs split in these two clusters
colors.single.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(single.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.single.X.scaled.data1)

# Have a look at the silhouette
sil.single.X.scaled.data1 <- silhouette(cutree(single.X.scaled.data1,2),dist.X.scaled.data1)
fviz_silhouette(sil.single.X.scaled.data1)

# Despite the ASW is big (0.64) this solution is awful because it only considers one cluster for all the points except for only one

##################################################################################################################
# Complete linkage
complete.X.scaled.data1<- hclust(dist.X.scaled.data1,method="complete")

# Plot dendogram of the solution and take k=4 as with Kmeans
plot(complete.X.scaled.data1,main="Complete linkage",cex=0.8)
rect.hclust(complete.X.scaled.data1,k=2,border="deepskyblue2")

# See the assignment
cutree(complete.X.scaled.data1,2)
table(cutree(complete.X.scaled.data1,2))

# Make a plot of the first two PCs split in these four clusters
colors.complete.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(complete.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.complete.X.scaled.data1)

# Have a look at the silhouette
sil.complete.X.scaled.data1 <- silhouette(cutree(complete.X.scaled.data1,2),dist.X.scaled.data1)
fviz_silhouette(sil.complete.X.scaled.data1)

# Better solution than the single linkage but it is still being not good since only few points are grouped in the second cluster and it has a lot of negative points

##################################################################################################################
# Average linkage
average.X.scaled.data1 <- hclust(dist.X.scaled.data1,method="average")

# Plot dendogram of the solution and take k=4 as with Kmeans
plot(average.X.scaled.data1,main="Average linkage",cex=0.8)
rect.hclust(average.X.scaled.data1,k=2,border="deepskyblue2")

# See the assignment
cutree(average.X.scaled.data1,2)
table(cutree(average.X.scaled.data1,2))

# Make a plot of the first two PCs split in these two clusters
colors.average.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(average.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.average.X.scaled.data1)

# Have a look at the silhouette
sil.average.X.scaled.data1 <- silhouette(cutree(average.X.scaled.data1,2),dist.X.scaled.data1)
fviz_silhouette(sil.average.X.scaled.data1)

# This solution is not good either, it happens the same than with single linkage

##################################################################################################################
# Centroids linkage
centroid.X.scaled.data1 <- hclust(dist.X.scaled.data1,method="centroid")

# Plot dendogram of the solution and take k=4 as with Kmeans
plot(centroid.X.scaled.data1,main="Centroid linkage",cex=0.8)
rect.hclust(centroid.X.scaled.data1,k=2,border="deepskyblue2")

# See the assignment
cutree(centroid.X.scaled.data1,2)
table(cutree(centroid.X.scaled.data1,2))

# Make a plot of the first two PCs split in these two clusters
colors.centroid.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(centroid.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.centroid.X.scaled.data1)

# Have a look at the silhouette
sil.centroid.X.scaled.data1 <- silhouette(cutree(centroid.X.scaled.data1,2),dist.X.scaled.data1)
fviz_silhouette(sil.centroid.X.scaled.data1)

# This solution is not good either, it happens the same than with single linkage

##################################################################################################################
# Ward method
ward.X.scaled.data1 <- hclust(dist.X.scaled.data1,method="ward.D2")

# Plot dendogram of the solution and take k=2 as with Kmeans
plot(ward.X.scaled.data1,main="Ward linkage",cex=0.8)
rect.hclust(ward.X.scaled.data1,k=2,border="deepskyblue2")

# See the assignment
cutree(ward.X.scaled.data1,2)
table(cutree(ward.X.scaled.data1,2))

# Make a plot of the first two PCs split in these two clusters
colors.ward.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(ward.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.ward.X.scaled.data1)

# Have a look at the silhouette
sil.ward.X.scaled.data1 <- silhouette(cutree(ward.X.scaled.data1,2),dist.X.scaled.data1)
fviz_silhouette(sil.ward.X.scaled.data1)

# This solution is slightly better, but it is not one of the highest (ASW=0.14)

##################################################################################################################
# Compute the Manhattan distance matrix between the observations in the data matrix

dist2.X.scaled.data1 <- daisy(X.scaled.data1,metric="manhattan",stand=FALSE)

##################################################################################################################
# Complete linkage
complete.X.scaled.data1<- hclust(dist2.X.scaled.data1,method="complete")

# Plot dendogram of the solution and take k=2 as with Kmeans
plot(complete.X.scaled.data1,main="Complete linkage",cex=0.8)
rect.hclust(complete.X.scaled.data1,k=2,border="deepskyblue2")

# See the assignment
cutree(complete.X.scaled.data1,2)
table(cutree(complete.X.scaled.data1,2))

# Make a plot of the first two PCs split in these four clusters
colors.complete.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(complete.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.complete.X.scaled.data1)

# Have a look at the silhouette
sil.complete.X.scaled.data1 <- silhouette(cutree(complete.X.scaled.data1,2),dist2.X.scaled.data1)
fviz_silhouette(sil.complete.X.scaled.data1)

# Better solution than with the euclidean distance but worse ASW (0.13) than the Kmeans 

##################################################################################################################
# Ward method
ward.X.scaled.data1 <- hclust(dist2.X.scaled.data1,method="ward.D2")

# Plot dendogram of the solution and take k=2 as with Kmeans
plot(ward.X.scaled.data1,main="Ward linkage",cex=0.8)
rect.hclust(ward.X.scaled.data1,k=2,border="deepskyblue2")

# See the assignment
cutree(ward.X.scaled.data1,2)
table(cutree(ward.X.scaled.data1,2))

# Make a plot of the first two PCs split in these two clusters
colors.ward.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(ward.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.ward.X.scaled.data1)

# Have a look at the silhouette
sil.ward.X.scaled.data1 <- silhouette(cutree(ward.X.scaled.data1,2),dist2.X.scaled.data1)
fviz_silhouette(sil.ward.X.scaled.data1)

# This solution is slightly better, even than the Kmeans (ASW=0.25), althought it still has a lot of negative points

##################################################################################################################
# Agglomerative hierarchical clustering analysis for the Census US data set - Gower (mixed data)!!!
##################################################################################################################

dist.Gower.scaled.data2

# Use ward method
complete.scaled.data2 <- hclust(dist.Gower.scaled.data2,method="ward.D2")

# Plot dendogram of the solution and take k=6
plot(complete.scaled.data2,main="Ward linkage",cex=0.8)
rect.hclust(complete.scaled.data2,k=6,border="deepskyblue2")

# See the assignment
cutree(complete.scaled.data2,6)
table(cutree(complete.scaled.data2,6))

# Make a plot of the first two PCs split in these six clusters
colors.complete.scaled.data2 <- c("deepskyblue2","firebrick2", "orange", "forestgreen", "black", "pink")[cutree(complete.scaled.data2,6)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.complete.scaled.data2)

# Have a look at the silhouette
sil.complete.scaled.data2 <- silhouette(cutree(complete.scaled.data2,6),dist.Gower.scaled.data2)
fviz_silhouette(sil.complete.scaled.data2)

summary(sil.complete.scaled.data2)

# ASW value of 0.08, it is the same than the one with Kmedoids (the partitionally method)

# Surely, a small number of clusters is more appropriate here, lets try with K=2 for comparing with Kmeans
################################################################

# Use Ward method
complete.scaled.data2 <- hclust(dist.Gower.scaled.data2,method="ward.D2")

# Plot dendogram of the solution and take k=2
plot(complete.scaled.data2,main="Ward linkage",cex=0.8)
rect.hclust(complete.scaled.data2,k=2,border="deepskyblue2")

# See the assignment
cutree(complete.scaled.data2,2)
table(cutree(complete.scaled.data2,2))

# Make a plot of the first two PCs split in these two clusters
colors.complete.scaled.data2 <- c("deepskyblue2","firebrick2")[cutree(complete.scaled.data2,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.complete.scaled.data2)

# Have a look at the silhouette
sil.complete.scaled.data2 <- silhouette(cutree(complete.scaled.data2,2),dist.Gower.scaled.data2)
fviz_silhouette(sil.complete.scaled.data2)

#The ASW is 0.09, it is better than k=6

##################################################################################################################
# Divisive hierarchical clustering analysis for the Census US data set
##################################################################################################################

Diana.X.scaled.data1 <- diana(X.scaled.data1,metric="manhattan")

# Plot dendogram of the solution
plot(Diana.X.scaled.data1,main="DIANA")

# Hit two times Return to see the dendrogram
# The heights here are the diameters of the clusters before splitting
# Take k=2

rect.hclust(Diana.X.scaled.data1,k=2,border="deepskyblue2")

# Make a plot of the first two PCs split in these two clusters
colors.Diana.X.scaled.data1 <- c("deepskyblue2","firebrick2")[cutree(Diana.X.scaled.data1,2)]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.Diana.X.scaled.data1)

# Have a look at the silhouette
sil.Diana.X.scaled.data1 <- silhouette(cutree(Diana.X.scaled.data1,2),dist.X.scaled.data1)
fviz_silhouette(sil.Diana.X.scaled.data1)

# The ASW value is 0.22, it is better than the mean of the values for others method althouth it is not the biggest, 
# nevertheless, in the graph we can see almost none negative points, so, this results looks like the best one

##################################################################################################################
##################################################################################################################
# Model-based clustering for the Census US data set
##################################################################################################################
##################################################################################################################

##################################################################################################################
# Install and load mclust package

install.packages("mclust")
library(mclust)

# Compute the value of the BIC for all the possible models

BIC.X.scaled.data1 <- mclustBIC(X.scaled.data1)
BIC.X.scaled.data1

# Note that the number of possible models is quite reduced because the number of observations is very small 
# compared with the dimension (n << p)

# Have a look at the different configurations

?mclustModelNames

# Have a look at the results

plot(BIC.X.scaled.data1)

# Note that the function returns -BIC, so we have to select the maximum value of this quantity

# Run Mclust for the optimal solution

X.scaled.data1.Mclust <- Mclust(X.scaled.data1,x=BIC.X.scaled.data1)

summary(X.scaled.data1.Mclust)

# Thus, Mclust selects a model with 7 clusters in which the covariance matrices are diagonal and have equal shapes

# The mixing probabilities can be found here
X.scaled.data1.Mclust$parameters$pro

# The sample means and covariance matrices can be found here
# Do not run in this case
X.scaled.data1.Mclust$parameters$mean
X.scaled.data1.Mclust$parameters$variance

# The probabilities of the observations belongs to each cluster
X.scaled.data1.Mclust$z

# The vector of clusters
X.scaled.data1.Mclust$classification

##################################################################################################################
# Have a look at the results with the first two PCs

colors.Mclust.X.scaled.data1 <- c("deepskyblue2","firebrick2","orange", "darkorchid1", "forestgreen", "black")[X.scaled.data1.Mclust$classification]
plot(PCS.X.scaled.data1$x[,1:2],pch=20,col=colors.Mclust.X.scaled.data1)

#Throught the plot it doesnt look good option, all the points are mixed



