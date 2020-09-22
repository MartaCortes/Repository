
rm(list = ls()) 

library(tidyverse)
library(MASS)
library(caret)
library(VGAM)
library(e1071) 
library(reshape) 
library(dplyr) 
library(ggplot2)
library(GGally)
library(psych)
library(gtable)
library(factoextra)
library(gridExtra)
library(fitdistrplus)
library(skimr)
library(mice)
library(VIM)
library(rpart)
library(pROC)
library(class)
library(randomForest)

# Loading and preparing data ----------------->>>>>>>>>>>>>

#The data classifies the leaves in 30 species, then, we have added type of "Deciduoud" or "Evergreen" of each specie in order 
# to have less groups, and predict that response. We do that because otherwise we will have very small probability 
# of belonging to a concrete group


data <- read_delim("~/Documents/Marta/MÁSTER DATA SCIENCE/Statistical Learning/Assignment Learning/VF/DataSet Leafs.csv", ";", escape_double = FALSE, trim_ws = TRUE)
data <- as.data.frame(data)

data[,17] <- NA
data[,17] <- ifelse(data[,1] == 1 | data[,1] == 4 | data[,1] == 8 | data[,1] == 16 |data[,1] == 17 |data[,1] == 18
                    |data[,1] == 19|data[,1] == 20|data[,1] == 21|data[,1] == 23|data[,1] == 24
                    |data[,1] == 25|data[,1] == 26|data[,1] == 28|data[,1] == 29|data[,1] == 30,"Evergreen",data[,17])
data[,17] <- ifelse(data[,1] == 2 | data[,1] == 3 | data[,1] == 5 | data[,1] == 6 |data[,1] == 7 |data[,1] == 9
                    |data[,1] == 10|data[,1] == 11|data[,1] == 12|data[,1] == 13|data[,1] == 14
                    |data[,1] == 15|data[,1] == 22|data[,1] == 27,"Deciduous",data[,17])


Leaves <- na.omit(data)
names(Leaves)[17] <- c("Category")

Leaves$Class=as.factor(Leaves$Class)
Leaves$`Specimen Number`=as.factor(Leaves$`Specimen Number`)
Leaves$Category=as.factor(Leaves$Category)



#### General Visualization ----------->>>>>>>>>>>>

# Bar chart about the distribution of the specimen numbers of each type
plot1 <- ggplot(Leaves, aes(x=Category,fill=(Category))) + geom_bar() + 
  xlab("Category") + ggtitle("Category distribution") +
  theme(legend.position="none")
plot1

# Densities and of each variable in order to analysis the distributions
Leaves <- Leaves[,-(1:2)]
Leaves_melt <- melt(Leaves, id="Category")

plot2 <- ggplot(Leaves_melt, aes (value, fill=variable)) +
  geom_density() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")
plot2

#Plots QQplot
plot3 <- ggplot(Leaves_melt, aes (sample=value, color=variable)) +
  geom_qq() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")
plot3


skim(Leaves)
summary(Leaves)
#Many of the variables might have some transformation
#Lets transform them --------------->>>>>>>>>>
library(rcompanion)
library(EnvStats)


Eccentricity2=(100-Leaves$Eccentricity)^2
qqPlot(Eccentricity2, add.line = TRUE)

`Aspect Ratio2` <- log(Leaves$`Aspect Ratio`+1)
qqPlot(`Aspect Ratio2`, add.line = TRUE)

Solidity2 <- (2-log(Leaves$Solidity))^-50
qqPlot(Solidity2, add.line = TRUE)

`Stochastic Convexity2` <- -1*(log(Leaves$`Stochastic Convexity`+1))^50
qqPlot(`Stochastic Convexity2`, add.line = TRUE)

`Isoperimetric Factor2` <- -1*(Leaves$`Isoperimetric Factor`)^1.5
qqPlot(`Isoperimetric Factor2`, add.line = TRUE)

`Maximal Indentation Depth2` <- log(Leaves$`Maximal Indentation Depth`)
qqPlot(`Maximal Indentation Depth2`, add.line = TRUE)

Lobeness2 <- log(8+log(Leaves$Lobedness))
qqPlot(Lobeness2, add.line = TRUE)

Smoothness2 <- log(Leaves$Smoothness)
qqPlot(Smoothness2, add.line = TRUE)

`Third moment2` <- log(Leaves$`Third moment`)
qqPlot(`Third moment2`, add.line = TRUE)

Uniformity2 <- -1*(1+log(Leaves$Uniformity))^4
qqPlot(Uniformity2, add.line = TRUE)

Entropy2 <- -1*(log(Leaves$Entropy-0.15))^2
qqPlot(Entropy2, add.line = TRUE)

Leaves2=mutate(Leaves, Eccentricity=(100-Eccentricity)^2, `Aspect Ratio`=log(`Aspect Ratio`+1), Solidity=(2-log(Solidity))^-50,
               `Stochastic Convexity`=-1*(log(`Stochastic Convexity`+1))^50, `Isoperimetric Factor`=-1*(`Isoperimetric Factor`)^1.5,
               `Maximal Indentation Depth`=log(`Maximal Indentation Depth`), Lobedness=log(8+log(Lobedness)), Smoothness=log(Smoothness),
               `Third moment`=log(`Third moment`), Uniformity=-1*(1+log(Uniformity))^4, Entropy=-1*(log(Entropy-0.15))^2)

#Checking again the distribution with the transformation done ------------>>>>>>>>>
Leaves2_melt <- melt(Leaves2)

plot4 <- ggplot(Leaves2_melt, aes (value, fill=variable)) +
  geom_density() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")
plot4

plot5 <- ggplot(Leaves2_melt, aes (sample=value, color=variable)) +
  geom_qq() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")
plot5

grid.arrange(plot2,plot4)
grid.arrange(plot3,plot5)
#Better

#Standarizations
scaled.Leaves2 <- as.data.frame(scale(Leaves2[,1:14]))

#Check that we get mean of 0 and sd of 1
colMeans(scaled.Leaves2) 
#Faster version of apply(scaled.dat, 2, mean)
apply(scaled.Leaves2, 2, sd)

scaled.Leaves2 <- cbind(Leaves2[,15],scaled.Leaves2)
names(scaled.Leaves2)[1]<-c("Category")

scaled.Leaves2_melt <- melt(scaled.Leaves2)

#Plots Boxplots
plot6<-ggplot(data = scaled.Leaves2_melt, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() +
  facet_wrap( ~ variable, scales="free") +
  theme(legend.position="none")
plot6


#Densities of each variable per each type of leaf
plot7 <- ggplot(scaled.Leaves2_melt, aes(x = value, fill=variable)) +  
  geom_density(aes(group = Category, 
                   colour = Category, 
                   fill = Category),
               alpha = 0.2) +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="bottom")
plot7


#Correlation Plot of all variables -------------------->>>>>>>>>>>>>

plot8 <- pairs.panels(scaled.Leaves2[2:15], 
                      method = "pearson", # correlation method
                      hist.col = "#00AFBB",
                      cex=5,
                      density = TRUE,  # show density plots
                      ellipses = TRUE # show correlation ellipses
)
plot8

#There are a lot of variables with high correlations, then, let's apply a dimensionality reduction with PCA, but first we have to scale the data.


####################
#PCA ##############
####################

#################################################################################################################
# Sample sizes and dimensions of the Leaves data set

X.LeavesPCA <- scaled.Leaves2[,2:15]

n.X.LeavesPCA <- nrow(X.LeavesPCA)
n.X.LeavesPCA
p.X.LeavesPCA <- ncol(X.LeavesPCA)
p.X.LeavesPCA

##################################################################################################################
# PCA for the Leaves data set
##################################################################################################################

##################################################################################################################
# There are several functions in R to perform PCA. Next, we use the function prcomp

PCS.X.LeavesPCA <- prcomp(X.LeavesPCA)

# Have a look at the outputs

names(X.LeavesPCA)

##################################################################################################################
# The PC scores is are in x. In this case, the scores are a matrix of size 64x64

dim(PCS.X.LeavesPCA$x)
head(PCS.X.LeavesPCA$x)

##################################################################################################################
# Make a plot of the first two PC scores
plot(PCS.X.LeavesPCA$x[,1:2],pch=19,col="deepskyblue2",main="First two PCs")
abline(h=0,v=0)

##################################################################################################################
# The eigenvalues of the sample covariance matrix of X, i.e., the variances of the PCs are the square of sdev
# As in this example n<p, only 14 eigenvalues can be different from 0. These are those that appear here.

PCS.X.LeavesPCA$sdev^2

##################################################################################################################
# Have a look at these eigenvalues

# Screeplot with the 14 eigenvalues

plot9<-fviz_eig(PCS.X.LeavesPCA,ncp=14,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")
plot9

# Screeplot with the first 5 eigenvalues
fviz_eig(PCS.X.LeavesPCA,ncp=5,addlabels=T,barfill="deepskyblue2",barcolor="deepskyblue4")

##################################################################################################################
# How many PCs are important?

# Have a look at the proportion of explained variance and the cumulative proportion of explained variance
get_eigenvalue(PCS.X.LeavesPCA)

# We only need 3 PCs to have more than 70% of the total variability of the data set explained

##################################################################################################################
# The loading matrix, i.e., the eigenvectors of the sample covariance matrix of X are given in rotation
# The output is a matrix of size 160x14, so we only have a look at the first few rows

dim(PCS.X.LeavesPCA$rotation)
head(PCS.X.LeavesPCA$rotation)

#Plot of first two PC scores with the initial variables
plot10<-biplot(PCS.X.LeavesPCA,col=c("deepskyblue2","firebrick2"),cex=c(0.5,0.8))
plot10

# This is useful to understand also which are the important variables in the data set in terms of variability
library(corrplot)
corrplot(cor(X.LeavesPCA,PCS.X.LeavesPCA$x),is.corr=T)

# Reduce to only three PCs
plot11<-corrplot(cor(X.LeavesPCA,PCS.X.LeavesPCA$x[,1:3]),is.corr=T)
plot11

PCS.X.LeavesPCA$rotation[,1:3]

#We can conclude that only 3 PCAs are enough for explain all the information of our data set and the varibles that contribute 
#the most to these PCAS are:
#Maximal Indentation Deep, Lobedness and Solidity for the first PCA
#Average Contrast, Smoothness and Entropy for the second PCA
#Eccentricity and Alongation for the third PCA
#Then, we are going to consider all the variables and these 8 variables instead of 14 for the classification modelling and compare both cases


###############################################################
#LOGISTIC REGRESSION

scaled.Leaves3<-scaled.Leaves2

names(scaled.Leaves3)[3]<-c("AspectRatio")
names(scaled.Leaves3)[6]<-c("StochasticConvexity")
names(scaled.Leaves3)[7]<-c("IsoperimetricFactor")
names(scaled.Leaves3)[8]<-c("MaximalIndentation")
names(scaled.Leaves3)[10]<-c("AverageIntensity")
names(scaled.Leaves3)[11]<-c("AverageContrast")
names(scaled.Leaves3)[13]<-c("ThirdMoment")

set.seed(0)
spl = createDataPartition(scaled.Leaves3$Category, p = 0.8, list = FALSE)  # 80% for training

LeavesTrain = scaled.Leaves3[spl,]
LeavesTest = scaled.Leaves3[-spl,]

str(LeavesTrain)
str(LeavesTest)

names(LeavesTrain)
dim(LeavesTrain)
summary(LeavesTrain)

#with ALL predictors
log.fit1 = vglm(Category ~., family=multinomial(refLevel=1), data=scaled.Leaves3)
summary(log.fit1)
# 2 groups, hence 1 regression

exp(coef(log.fit1))

# Interpretation:
#The variables which the exponential of coefficients (odds) is higher than 1, represent better the Deciduous leaves
#since the reference level is Deciduous. Then, the odds which are less than one represent better the Evergreen leaves.


# predicting the testing set
# output are probabilities, no labels
prob.test1 <- predict(log.fit1,newdata = LeavesTest,type = "response")
head(prob.test1)

# We can apply Bayes rule: maximum probability
pred.test1 <- ifelse(prob.test1 > 0.5,1,0)
head(pred.test1)

# summarize accuracy (confusion matrix) for a given probability rule (maximum in this case)
# predictions in rows, true values in columns (but we can change the order)
ConfMat1 = table(pred.test1,LeavesTest$Category)
ConfMat1

n = length(LeavesTest$Category)
prop.errors1 <- (n - sum(diag(ConfMat1))) / n
prop.errors1
#Error of 25%

accuracy1 <- sum(diag(ConfMat1)) / n
accuracy1
#Accuracy of 75%
# Reasonable error and accuracy

# With the EIGTH predictors of the PCA:
log.fit2 = vglm(Category ~ Alongation+Lobedness+Solidity+MaximalIndentation+AverageContrast+Smoothness+Entropy+Eccentricity, family=multinomial(refLevel=1), data=LeavesTrain)
summary(log.fit2)

prob.test2 = predict(log.fit2, newdata=LeavesTest, type="response")
head(prob.test2)
# output are probabilities, no labels

# We can apply Bayes rule: maximum probability
pred.test2 <- ifelse(prob.test2 > 0.5,1,0)
head(pred.test2)

# summarize accuracy (confusion matrix) for a given probability rule (maximum in this case)
# predictions in rows, true values in columns (but we can change the order)
ConfMat2 = table(pred.test2,LeavesTest$Category)
ConfMat2

n = length(LeavesTest$Category)
prop.errors2 <- (n - sum(diag(ConfMat2))) / n
prop.errors2
#Error of 36%

accuracy2 <- sum(diag(ConfMat2)) / n
accuracy2
#Accuracy of 64%

# Worse results with the 8 predictors of the PCA than with all the variables (completed model)

#Lets try with some interaction
log.fit3 = vglm(Category ~ Alongation+Lobedness+Solidity+MaximalIndentation+AverageContrast+Smoothness+Entropy+Eccentricity+AspectRatio+StochasticConvexity+IsoperimetricFactor+AverageIntensity+ThirdMoment+Uniformity+Solidity:AverageIntensity, family=multinomial(refLevel=1), data=LeavesTrain)
summary(log.fit3)

prob.test3 = predict(log.fit3, newdata=LeavesTest, type="response")
head(prob.test3)
# output are probabilities, no labels

# We can apply Bayes rule: maximum probability
pred.test3 <- ifelse(prob.test3 > 0.5,1,0)
head(pred.test3)

# summarize accuracy (confusion matrix) for a given probability rule (maximum in this case)
# predictions in rows, true values in columns (but we can change the order)
ConfMat3 = table(pred.test3,LeavesTest$Category)
ConfMat3

n = length(LeavesTest$Category)
prop.errors3 <- (n - sum(diag(ConfMat3))) / n
prop.errors3
#Error of 32%

accuracy3 <- sum(diag(ConfMat3)) / n
accuracy3
#Accuracy of 67%
#No better with interactions
#Best model: log.fit1


########################################
################### BAYES CLASSIFICATION

#Although we already tried to transform the data to achieve normality, we didn't get it at all, therefore, we need
#to assume normality for doing LDA and QDA, even if we don't have it at all

#Matrix of variances  --------------------->>>>>>>>>>>>>>
n <- dim(scaled.Leaves3[,2:15])[1]
n
p <- dim(scaled.Leaves3[,2:15])[2]
p


LeavesBayesVar <- group_by(scaled.Leaves3,Category) %>% summarise_all(var)
LeavesBayesVar

#We do not have same variances sometimes, although they are not too different
#In any case, we also need to assume equal variances for doing LDA and QDA
#The only consideration is that we will not have as good prediction results as possible because we have errors 
#in the assumptions (normality and equal variances)


#With ALL variables
# LDA
lda.model1 <- lda(Category ~ ., data=LeavesTrain)
lda.model1

probability1 = predict(lda.model1, newdata=LeavesTest)$posterior
head(probability1)

# We apply the Bayes rule of maximum probability
prediction1 <- max.col(probability1)
head(prediction1)

# It's equivalent to
prediction1 = predict(lda.model1, newdata=LeavesTest)$class
head(prediction1)

# Performance: confusion matrix
confusionMatrix(prediction1, LeavesTest$Category)

ConfMat4 = table(prediction1,LeavesTest$Category)
ConfMat4

n = length(LeavesTest$Category)
prop.errors4 <- (n - sum(diag(ConfMat4))) / n
prop.errors4
#Error of 27%

accuracy4 <- sum(diag(ConfMat4)) / n
accuracy4
#Accuracy of 73%


# Quadratic Discriminant Analysis (QDA)
qda.model1 <- qda(Category ~ ., data=LeavesTrain)
qda.model1

prediction2 = predict(qda.model1, newdata=LeavesTest)$class
confusionMatrix(prediction2, LeavesTest$Category)

ConfMat5 = table(prediction2,LeavesTest$Category)
ConfMat5

n = length(LeavesTest$Category)
prop.errors5 <- (n - sum(diag(ConfMat5))) / n
prop.errors5
#Error of 22%

accuracy5 <- sum(diag(ConfMat5)) / n
accuracy5
#Accuracy of 78%


# Naive Bayes (Gaussian and linear)
naive.model1 <- naiveBayes(Category ~ ., data=LeavesTrain, laplace=1)
# Laplace smoothing is used to avoid 0-class probabilities
naive.model1

prediction3 = predict(naive.model1, newdata=LeavesTest)
confusionMatrix(prediction3, LeavesTest$Category)

ConfMat6 = table(prediction3,LeavesTest$Category)
ConfMat6

n = length(LeavesTest$Category)
prop.errors6 <- (n - sum(diag(ConfMat6))) / n
prop.errors6
#Error of 27%

accuracy6 <- sum(diag(ConfMat6)) / n
accuracy6
#Accuracy of 73%


#With the 8 selected variables of the PCA
# LDA
lda.model2 <- lda(Category ~ Alongation+Lobedness+Solidity+MaximalIndentation+AverageContrast+Smoothness+Entropy+Eccentricity, data=LeavesTrain)
lda.model2

probability2 = predict(lda.model2, newdata=LeavesTest)$posterior
head(probability2)

# We apply the Bayes rule of maximum probability
prediction4 <- max.col(probability2)
head(prediction4)

# It's equivalent to
prediction4 = predict(lda.model2, newdata=LeavesTest)$class
head(prediction4)

# Performance: confusion matrix
confusionMatrix(prediction4, LeavesTest$Category)

ConfMat7 = table(prediction4,LeavesTest$Category)
ConfMat7

n = length(LeavesTest$Category)
prop.errors7 <- (n - sum(diag(ConfMat7))) / n
prop.errors7
#Error of 28%

accuracy7 <- sum(diag(ConfMat7)) / n
accuracy7
#Accuracy of 72%


# Quadratic Discriminant Analysis (QDA)
qda.model2 <- qda(Category ~ Alongation+Lobedness+Solidity+MaximalIndentation+AverageContrast+Smoothness+Entropy+Eccentricity, data=LeavesTrain)
qda.model2

prediction5 = predict(qda.model2, newdata=LeavesTest)$class
confusionMatrix(prediction5, LeavesTest$Category)

ConfMat8 = table(prediction5,LeavesTest$Category)
ConfMat8

n = length(LeavesTest$Category)
prop.errors8 <- (n - sum(diag(ConfMat8))) / n
prop.errors8
#Error of 24%

accuracy8 <- sum(diag(ConfMat8)) / n
accuracy8
#Accuracy of 76%


# Naive Bayes (Gaussian and linear)
naive.model2 <- naiveBayes(Category ~ Alongation+Lobedness+Solidity+MaximalIndentation+AverageContrast+Smoothness+Entropy+Eccentricity, data=LeavesTrain, laplace=1)
# Laplace smoothing is used to avoid 0-class probabilities
naive.model2

prediction6 = predict(naive.model2, newdata=LeavesTest)
confusionMatrix(prediction6, LeavesTest$Category)

ConfMat9 = table(prediction6,LeavesTest$Category)
ConfMat9

n = length(LeavesTest$Category)
prop.errors9 <- (n - sum(diag(ConfMat9))) / n
prop.errors9
#Error of 25%

accuracy9 <- sum(diag(ConfMat9)) / n
accuracy9
#Accuracy of 75%

#The best models are the ones with all predictors and then the QDA method is the best, next is the Logistic Regression 
#and then, the Naive Bayes and the last one is the LDA


#############################################
############# MACHINE-LEARNING

###############################
################# KNN

knn_pred1 <- knn(train=LeavesTrain[,-1], test=LeavesTest[,-1], cl = LeavesTrain$Category, k=3)
confusionMatrix(knn_pred1, LeavesTest$Category)
ConfMat10 = table(knn_pred1,LeavesTest$Category)
n = length(LeavesTest$Category)
prop.errors10 <- (n - sum(diag(ConfMat10))) / n
prop.errors10
#Accuracy: 0.69, Kappa: 0.38, error: 0.31

#In essence, the kappa statistic is a measure of how closely the instances classified by the machine 
#learning classifier matched the data labeled as ground truth, controlling for the accuracy of a random 
#classifier as measured by the expected accuracy.

# how to choose the hyper-parameter k? 

ctrl <- trainControl(method = "cv", number = 5,
                     classProbs = TRUE,
                     verboseIter=T)

result_knnFit_Accuracy<-c(0)
result_knnFit_Kappa<-c(0)
for (i in 1:30) {
knnFit <- train(Category ~ ., 
                method = "knn", 
                data = LeavesTrain,
                tuneLength = 3,
                metric = "Accuracy",
                trControl = ctrl)
print(knnFit)
knnPred = predict(knnFit, LeavesTest)
result_knnFit_Accuracy[i]<-confusionMatrix(knnPred,LeavesTest$Category)$overall[1]
result_knnFit_Kappa[i]<-confusionMatrix(knnPred,LeavesTest$Category)$overall[2]
}
mean(result_knnFit_Accuracy)
mean(result_knnFit_Kappa)
prop.errors11 <- (1-result_knnFit_Accuracy)
mean(prop.errors11)
#Due to the random initialization of the weights the final results can differ a lot, so, we have take the most repeated value
#of k, and we have calculated the result with it. The optimal k value is k=7. Then, the outputs are accuracy: 0.78, Kappa:0.56,
#error:0.22

#A bit better with the tunning

###############################
##################### SVM

svm.train <- svm(Category ~., data=LeavesTrain, scale=F, kernel="radial",
                 gamma=0.01, cost=1)
# gamma is the parameter of the radial basis function (kernel)
svm.pred1 <- predict(svm.train, newdata=LeavesTest)
confusionMatrix(svm.pred1, LeavesTest$Category)
ConfMat11 = table(svm.pred1, LeavesTest$Category)
n = length(LeavesTest$Category)
prop.errors12 <- (n - sum(diag(ConfMat11))) / n
prop.errors12
#Accuracy: 0.75, Kappa:0.49, error:0.25

#how to select the hyper-parameters gamma and cost?

result_svmFit_Accuracy<-c(0)
result_svmFit_Kappa<-c(0)
for (i in 1:30) {
svmFit <- train(Category ~., method = "svmRadial", 
                data = LeavesTrain,
                tuneGrid = expand.grid(C = c(.25, .5, 1),
                                       sigma = c(0.01,.05)), 
                metric = "Accuracy",
                trControl = ctrl)
print(svmFit)
svmPred = predict(svmFit, LeavesTest)
result_svmFit_Accuracy[i]<-confusionMatrix(svmPred,LeavesTest$Category)$overall[1]
result_svmFit_Kappa[i]<-confusionMatrix(svmPred,LeavesTest$Category)$overall[2]
}
mean(result_svmFit_Accuracy)
mean(result_svmFit_Kappa)
prop.errors13 <- (1-result_svmFit_Accuracy)
mean(prop.errors13)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy: 0.81, Kappa:0.61, error:0.19

#Much better with the hyper-parameter tunning and better than the KNN

##################################
############### Decision Trees

rpart.fit1 <- rpart(Category ~., method="class", data = LeavesTrain)
summary(rpart.fit1)
# For each node in the tree, the number of examples reaching the decision point is listed

# Visualizing decision trees
library(rpart.plot)
plot11<-rpart.plot(rpart.fit1, digits = 3, fallen.leaves = TRUE,
           type = 3, extra=101)
plot11

DTpred1 <- predict(rpart.fit1, LeavesTest, type="class")
confusionMatrix(DTpred1, LeavesTest$Category)
ConfMat12 = table(svm.pred1, LeavesTest$Category)
n = length(LeavesTest$Category)
prop.errors14 <- (n - sum(diag(ConfMat12))) / n
prop.errors14
#Accuracy: 0.77, Kappa:0.55, error:0.25

# C5.0: advanced DT model, somehow it's a boosting approach 

result_C5.0_Accuracy<-c(0)
result_C5.0_Kappa<-c(0)
for (i in 1:30) {
fit.c50 <- train(Category ~.,
                 data=LeavesTrain,
                 method="C5.0",
                 metric="Accuracy",
                 tuneGrid = expand.grid(.winnow = c(TRUE,FALSE),.trials=c(15),.model="tree"),
                 trControl = ctrl)
print(fit.c50)
# Winnowing is a feature selection step conducted before modeling
# trials = number of boosting iterations
c50.pred <- predict(fit.c50, newdata=LeavesTest)
result_C5.0_Accuracy[i]<-confusionMatrix(c50.pred, LeavesTest$Category)$overall[1]
result_C5.0_Kappa[i]<-confusionMatrix(c50.pred, LeavesTest$Category)$overall[2]
}
mean(result_C5.0_Accuracy)
mean(result_C5.0_Kappa)
prop.errors15 <- (1-result_C5.0_Accuracy)
mean(prop.errors15)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy:0.78, kappa:0.56 and error:0.22

#Same result with tunning but better than KNN. Worse than SVM

###############################
############### RF

result_rf1_Accuracy<-c(0)
result_rf1_Kappa<-c(0)
for (i in 1:30) {
rf.train1 <- randomForest(Category ~., data=LeavesTrain,
                         ntree=250,mtry=4,cutoff=c(0.5,0.5),importance=TRUE, do.trace=T)
rf.pred1 <- predict(rf.train1, newdata=LeavesTest)
result_rf1_Accuracy[i]<-confusionMatrix(rf.pred1, LeavesTest$Category)$overall[1]
result_rf1_Kappa[i]<-confusionMatrix(rf.pred1, LeavesTest$Category)$overall[2]
}
mean(result_rf1_Accuracy)
mean(result_rf1_Kappa)
prop.errors16 <- (1-result_rf1_Accuracy)
mean(prop.errors16)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy:0.75, kappa:0.51 and error:0.25

# the cutoff in RF controls the probability to belong a group, since there is only two categories, it should be (0.5,0.5)
# but it works in the same way as the threshold probability in Bayes rule

result_rf2_Accuracy<-c(0)
result_rf2_Kappa<-c(0)
for (i in 1:30) {
rf.train2 <- train(Category ~., 
                  method = "rf", 
                  data = LeavesTrain,
                  ntree = 250,
                  cutoff=c(0.5,0.5),
                  tuneGrid = expand.grid(mtry = c(2,4,6,8,10,12,14)), 
                  metric = "Accuracy",
                  trControl = ctrl)
print(rf.train2)
rf.pred <- predict(rf.train2, newdata=LeavesTest)
result_rf2_Accuracy[i]<-confusionMatrix(rf.pred1, LeavesTest$Category)$overall[1]
result_rf2_Kappa[i]<-confusionMatrix(rf.pred1, LeavesTest$Category)$overall[2]
}
mean(result_rf2_Accuracy)
mean(result_rf2_Kappa)
prop.errors17 <- (1-result_rf2_Accuracy)
mean(prop.errors17)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy: 0.81, kappa:0.61 and error:0.19

#A bit better with the tunning, KNN and DT, and worse than SVM

####################################
####### Basic ensemble prediction

# Create mode function
mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

lda.model <- lda(Category ~ ., data=LeavesTrain)
ldaProb = predict(lda.model, newdata=LeavesTest)$posterior
threshold = 0.5
lda.pred = rep("Deciduous", nrow(LeavesTest))
lda.pred[which(ldaProb[,2] > threshold)] = "Evergreen"
confusionMatrix(factor(lda.pred), LeavesTest$Category)
ConfMat13 = table(lda.pred, LeavesTest$Category)
n = length(LeavesTest$Category)
prop.errors18 <- (n - sum(diag(ConfMat13))) / n
prop.errors18
#Accuracy: 0.73, Kappa:0.46, error:0.27

#Not one of the best results


result_ensemble_Accuracy<-c(0)
result_ensemble_Kappa<-c(0)
for (i in 1:30) {
ensemble.pred = apply(data.frame(lda.pred, knnPred, svmPred, c50.pred, rf.pred), 1, mode) 
result_ensemble_Accuracy[i]<-confusionMatrix(factor(ensemble.pred), LeavesTest$Category)$overall[1]
result_ensemble_Kappa[i]<-confusionMatrix(factor(ensemble.pred), LeavesTest$Category)$overall[2]
}
mean(result_ensemble_Accuracy)
mean(result_ensemble_Kappa)
prop.errors19 <- (1-result_ensemble_Accuracy)
mean(prop.errors19)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy:0.79, kappa:0.58 and error:0.21

#Not better than SVM at all


#### OTHER MODELS

###############################
######### Gradient Boosting
library(gbm)

result_GBM_Accuracy<-c(0)
result_GBM_Kappa<-c(0)
for (i in 1:30) {
GBM.train <- gbm(ifelse(LeavesTrain$Category=="Deciduous",0,1) ~., data=LeavesTrain,
                 distribution= "bernoulli",n.trees=250,shrinkage = 0.01,interaction.depth=2,n.minobsinnode = 8)
threshold = 0.5
gbmProb = predict(GBM.train, newdata=LeavesTest, n.trees=250, type="response")
gbmPred = rep("Deciduous", nrow(LeavesTest))
gbmPred[which(gbmProb > threshold)] = "Evergreen"
result_GBM_Accuracy[i]<-confusionMatrix(factor(gbmPred), LeavesTest$Category)$overall[1]
result_GBM_Kappa[i]<-confusionMatrix(factor(gbmPred), LeavesTest$Category)$overall[2]
}
mean(result_GBM_Accuracy)
mean(result_GBM_Kappa)
prop.errors20 <- (1-result_GBM_Accuracy)
mean(prop.errors20)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy:0.71, kappa:0.42 and error:0.29

#Quite bad...


# Now optimize the hyper-parameters with Caret:

xgb_grid = expand.grid(
  nrounds = c(500,1000),
  eta = c(0.01, 0.001), # c(0.01,0.05,0.1)
  max_depth = c(2, 4, 6),
  gamma = 1,
  colsample_bytree = c(0.2, 0.4),
  min_child_weight = c(1,5),
  subsample = 1
)

result_xgb_Accuracy<-c(0)
result_xgb_Kappa<-c(0)
for (i in 1:30) {
xgb.train = train(Category ~ .,  data=LeavesTrain,
                  trControl = ctrl,
                  metric="Accuracy",
                  maximize = F,
                  tuneGrid = xgb_grid,
                  method = "xgbTree"
)

# Variable importance
xgb_imp <- varImp(xgb.train, scale = F)
plot12<-plot(xgb_imp, scales = list(y = list(cex = .95)))
threshold = 0.5
xgbProb = predict(xgb.train, newdata=LeavesTest, type="prob")
xgbPred = rep("Deciduous", nrow(LeavesTest))
xgbPred[which(xgbProb[,2] > threshold)] = "Evergreen"
result_xgb_Accuracy[i]<-confusionMatrix(factor(xgbPred), LeavesTest$Category)$overall[1]
result_xgb_Kappa[i]<-confusionMatrix(factor(xgbPred), LeavesTest$Category)$overall[2]
}
plot12
mean(result_xgb_Accuracy)
mean(result_xgb_Kappa)
prop.errors21 <- (1-result_xgb_Accuracy)
mean(prop.errors21)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy:0.77, kappa:0.55 and error:0.23

#Better result with the tunning but worse than previous models


###############################
######### Neural Networks

ctrl$sampling <- NULL

result_nn_Accuracy<-c(0)
result_nn_Kappa<-c(0)
for (i in 1:30) {
nn.train <- train(Category ~., 
                  method = "avNNet",
                  repeats=15,
                  data = LeavesTest,
                  maxit = 1000,
                  tuneGrid = expand.grid(size=c(2,4,6), decay=c(0.01,0.001), bag=F), 
                  metric = "Accuracy",
                  trControl = ctrl)
plot13<-plot(nn.train)
nn_imp <- varImp(nn.train, scale = F)
plot14<-plot(nn_imp, scales = list(y = list(cex = .95)))
threshold = 0.5
nnProb = predict(nn.train, newdata=LeavesTest, type="prob")
nnPred = rep("Deciduous", nrow(LeavesTest))
nnPred[which(nnProb[,2] > threshold)] = "Evergreen"
result_nn_Accuracy[i]<-confusionMatrix(factor(nnPred), LeavesTest$Category)$overall[1]
result_nn_Kappa[i]<-confusionMatrix(factor(nnPred), LeavesTest$Category)$overall[2]
}
plot13
plot12
mean(result_nn_Accuracy)
boxplot(result_nn_Accuracy)
mean(result_nn_Kappa)
prop.errors22 <- (1-result_nn_Accuracy)
mean(prop.errors22)
#Due to the random initialization of the weights the final results can differ a lot, but as we can see in the result
#vector, the means of the outputs are accuracy:0.97, kappa:0.95 and error:0.03
 


