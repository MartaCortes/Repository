# Script created to do PRS
#       
# ==============================================
# Description - POLYGENIC RISK CALCULATION
# ==============================================


### INTELLIGENCE

library(base)
setwd("~/Documents/Marta/MASTER_DATA_SCIENCE/TFM/Genéticos/Código/Inteligencia")

# Read summary statistics
dataDiscovery <- read.table('SavageJansen_2018_intelligence_metaanalysis.txt', header=TRUE, stringsAsFactors = FALSE)
head(dataDiscovery)


### A) FILTRO DISCOVERY ###

#filtro de variantes imputadas a mala calidad
temp_discovery1 <- subset(dataDiscovery, minINFO>=0.8)   # SI EXISTE

#Eliminamos covariantes con strand ambiguity. Hay que comprobar qué columna es cromosoma, posición genómica, y cuales alelo 1 y 2.
#names(temp_discovery1)[4] <- c('BP')
temp_discovery1$POS <- paste(temp_discovery1$CHR, temp_discovery1$POS, sep=":")
temp_discovery1$C <- paste(temp_discovery1$A2, temp_discovery1$A1, sep="")
temp_discovery2 <- subset(temp_discovery1, !(C %in% c("at","ta", "cg", "gc")))
#temp_discovery2 <- temp_discovery2[,-2]


#Eliminamos también indels o trialélicas. 
temp_discovery3 <- subset(temp_discovery2, (C %in% c("ac", "ag", "tg", "tc", "ga", "gt", "ca", "ct")))
temp_discovery3$A1 = toupper(temp_discovery3$A1)
temp_discovery3$A2 = toupper(temp_discovery3$A2)
temp_discovery3$C = toupper(temp_discovery3$C)


#creo un archivo con frecuencia de cada variante por si aparece repetida, y eliminamos repetidas
FRECUENCIAS_discovery <- as.data.frame(table(temp_discovery3$POS))
temp_discovery3_UNIQ <- temp_discovery3[temp_discovery3$POS %in% FRECUENCIAS_discovery$Var1[FRECUENCIAS_discovery$Freq < 2], ]
dataDiscovery <- temp_discovery3_UNIQ
head(dataDiscovery)


### B) FILTRO TARGET ##

## USO SOLO VARIANTES FUERA DEL MHC. Leo ya el archivo de PGC modificado para ello.
setwd("~/Documents/Marta/MASTER_DATA_SCIENCE/TFM/Genéticos/Código/Inteligencia")
dataTarget <- read.table('EUGEI_WP2_WP6_PassingQC_AllPopulations.bim', header = FALSE, stringsAsFactors = FALSE)


# Change column names
colnames(dataTarget) <- c("V1", "V2", "V3", "POS", "A1", "A2")
dataTarget

# Creo una columna con la variable toto chr:pos; y ordeno
dataTarget$varID <- paste(dataTarget$V1, dataTarget$POS, sep=":")
dataTarget_1 <- dataTarget[,c("V1", "varID", "V3", "POS", "A1", "A2")]
dataTarget<-dataTarget_1

# Exportamos el archivo dataTarget y los sustituímos por el .bim de EUGEI (para que tenga la columna de variantes como CHR:POS)
write.table(dataTarget, 'dataTarget.txt', sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#creo un archivo con frecuencia de cada variante por si aparece repetida, y eliminamos repetidas
FRECUENCIAS_target<-as.data.frame(table(dataTarget$varID))
temp_target_UNIQ<-dataTarget[dataTarget$varID %in% FRECUENCIAS_target$Var1[FRECUENCIAS_target$Freq < 2], ]
dataTarget <- temp_target_UNIQ
head(dataTarget)
head(dataDiscovery)


### hacemos OVERLAP
# Filter the dataTarget by the variants in dataDiscovery
dataTarget_filtered <- dataTarget[which(dataTarget$varID %in% dataDiscovery$POS),]
# Filter the dataDiscovery by the variants in dataTarget_filtered
dataDiscovery_filtered <- dataDiscovery[which(dataDiscovery$POS %in% dataTarget_filtered$varID),]

rm(temp_discovery1)
rm(temp_discovery2)
rm(temp_discovery3)
rm(temp_discovery3_UNIQ)


### C) CLUNPING - ELIMINAR VARIANTES EN LD

setwd("~/Documents/Marta/MASTER_DATA_SCIENCE/TFM/Genéticos/Código/Inteligencia")

# preparar los archivos de discovery, que se utilizarán la P para el clump
# Me quedo con las columnas que me interesan
dataDiscovery_filtered1<-dataDiscovery_filtered[,c("CHR", "POS", "A1", "A2", "EAF_HRC", "stdBeta", "SE", "P")]
# cambio los nombres de encabezado
colnames(dataDiscovery_filtered1)<-c("CHR", "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P")
dataDiscovery_filtered <- dataDiscovery_filtered1

#head(dataDiscovery_filtered1)
write.table(dataDiscovery_filtered, "discovery_paraCLUMP.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Get the common variants and save a txt
varIDs <- data.frame(dataTarget_filtered$varID)
# write.table(dataDiscovery_filtered, 'dataDiscovery_filtered.txt', row.names = FALSE, col.names = TRUE)
write.table(varIDs, 'variableIDs_forFiltering.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)

# Create new .bed .bim and .fam files with only the variants in common among discovery and target
system2("./plink", args=c("--bfile EUGEI_WP2_WP6_PassingQC_AllPopulations", "--extract variableIDs_forFiltering.txt", "--make-bed", "--out EUGEI_WP2_WP6_new_filtered"))

# Hacemos CLUMP en plink con el archivo modificado de discovery
system2("./plink", args=c("--bfile EUGEI_WP2_WP6_new_filtered", "--clump discovery_paraCLUMP.txt", "--clump-p1 1", "--clump-p2 1", "--clump-kb 500", "--clump-r2 0.1", "--out EUGEI_WP2_WP6_new_filtered_CLUMPED"))
CLUMPED_PLINK<-read.table("EUGEI_WP2_WP6_new_filtered_CLUMPED.clumped", header=T)
write.table(CLUMPED_PLINK$SNP, "CLUMPED_PLINK.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
system2("./plink", args=c("--bfile ~/Documents/Marta/MASTER_DATA_SCIENCE/TFM/Genéticos/Código/Inteligencia/EUGEI_WP2_WP6_new_filtered", "--extract CLUMPED_PLINK.txt", "--make-bed", "--out EUGEI_WP2_WP6_new_filtered_new_CLUMPED"))

#rm(temp*)




### D)  CALCULO DE PRS SCORES  #########
head(dataDiscovery_filtered)
head(dataTarget_filtered)


#CREAMOS LOS ARCHIVOS NECESARIOS
SCZ_SCORE <- dataDiscovery_filtered[,c("SNP", "A1", "BETA")]
SCZ_P <- dataDiscovery_filtered[,c("SNP", "P")]
write.table(SCZ_SCORE, "SCZ_SCORE.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(SCZ_P, "SCZ_P.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Calculate PRS
system2("./plink", args=c("--bfile ~/Documents/Marta/MASTER_DATA_SCIENCE/TFM/Genéticos/Código/Inteligencia/EUGEI_WP2_WP6_new_filtered_new_CLUMPED", "--score SCZ_SCORE.txt", "--q-score-file SCZ_P.txt", "--q-score-range Q_RANGES.txt", "--out Intelligence_PRS_SCORE"))


#### PCA calculation
system2("./plink", args=c("--bfile EUGEI_WP2_WP6_PassingQC_AllPopulations", "--indep-pairwise 500 1 0.1", "--out PCA_PRUNED"))
system2("./plink", args=c("--bfile EUGEI_WP2_WP6_PassingQC_AllPopulations", "--extract PCA_PRUNED.prune.in", "--make-bed", "--out EUGEI_WP2_WP6_PRUNED"))
system2("./plink", args=c("--bfile EUGEI_WP2_WP6_PRUNED", "--pca", "--out PCA"))


#########
## REGRESSION GWAS MONTHLY INCOME
########


rm(list=ls())


library(readxl)
library(xlsx)
library(MASS)
library(dplyr)
library(ggplot2)
library(GGally)
library(psych)
library(reshape2)
library(ISLR)
library(factoextra)
library(plotrix)
library(gridExtra)
library(grid)
library(effects)
library(carData)
library(Epi)
library(pROC)
library(caret)
library(ordinalForest)
library(ggplot2)
library(glmnet)
library(skimr)
library(cbar)


setwd("~/Documents/Marta/MASTER_DATA_SCIENCE/TFM/Genéticos/Código")
table <- read_excel("Final Regression Dataset.xlsx")


###########################################  PREPROCESSING -----------------------------------------------------------

table <- table[,-c(1:2)]
table <- table[table[,1]==3,]
table <- table[,-1]


names(table)[4] <- c("EduLevel")
names(table)[5] <- c("SocCla")
names(table)[6] <- c("Income")
names(table)[7] <- c("ChildTrauma")
names(table)[8] <- c("Cannabis")
names(table)[19] <- c("ScoreInt")
names(table)[20] <- c("ScoreCog")
names(table)[21] <- c("ScoreEA")

table[table=="NA"] <- NA

table[,1] <- as.factor(unlist((table[,1])))
table[,2] <- as.numeric(unlist(table[,2]))
table[,3] <- as.factor(unlist(table[,3]))
table[,4] <- as.factor(unlist(table[,4]))
table[,5] <- as.factor(unlist(table[,5]))
table[,6] <- as.numeric(unlist(table[,6]))
table[,7] <- as.factor(unlist(table[,7]))
table[,8] <- as.factor(unlist(table[,8]))


for (i in 9:21) {
  table[,i] <- as.numeric(unlist(table[,i]))
}

#Grouping the Phenotypes in only 3 groups

table[,4] <- as.numeric(unlist(table[,4]))

table$EduLevel <- ifelse(table$EduLevel == 1 | table$EduLevel == 2, '1', table$EduLevel)
table$EduLevel <- ifelse(table$EduLevel == 3 | table$EduLevel == 4, '2', table$EduLevel)
table$EduLevel <- ifelse(table$EduLevel == 5 | table$EduLevel == 6, '3', table$EduLevel)

table[,4] <- as.factor(unlist(table[,4]))


table[,5] <- as.numeric(unlist(table[,5]))

table$SocCla <- ifelse(table$SocCla == 1 | table$SocCla == 2, '1', table$SocCla)
table$SocCla <- ifelse(table$SocCla == 3 | table$SocCla == 4, '2', table$SocCla)
table$SocCla <- ifelse(table$SocCla == 5 | table$SocCla == 6, '3', table$SocCla)

table[,5] <- as.factor(unlist(table[,5]))


#Visualization of all numeric variables (before deleting missing values)
table_melt <- melt(table[,-c(1,3:5,7:8)])

#Density Plots
p7 <- ggplot(table_melt, aes (value, fill=variable)) +
  geom_density() +
  facet_wrap(~variable, scales="free") +
  theme(legend.position="none")
#Maybe I have to apply a transformation in the Income variable


#Boxplots
p8 <- ggplot(table_melt, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() +
  facet_wrap( ~ variable, scales="free") +
  theme(legend.position="none")

grid.arrange(p7, p8, ncol = 2)

#Categorical variables distribution

Sex <- table[,1]
sum(is.na(Sex))
Sex <- na.omit(Sex)
Sex <- Sex %>% dplyr::count(Sex)

p1 <- ggplot(Sex, aes(x = Sex, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") + 
  xlab("Sex")+
  ylab("Count")+
  theme(axis.text.x = element_text(size=7,angle=0, color="skyblue4"))


Country <- table[,3]
sum(is.na(Country))
Country <- na.omit(Country)
Country <- Country %>% dplyr::count(Country)

p2 <- ggplot(Country, aes(x = Country, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") + 
  xlab("Country")+
  ylab("Count")+
  theme(axis.text.x = element_text(size=7,angle=0, color="skyblue4"))


EduLevel <- table[,4]
sum(is.na(EduLevel))
EduLevel <- na.omit(EduLevel)
EduLevel <- EduLevel %>% dplyr::count(EduLevel)

p3 <- ggplot(EduLevel, aes(x = EduLevel, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") + 
  xlab("EduLevel")+
  ylab("Count")+
  theme(axis.text.x = element_text(size=7,angle=0, color="skyblue4"))


SocCla <- table[,5]
sum(is.na(SocCla))
SocCla <- na.omit(SocCla)
SocCla <- SocCla %>% dplyr::count(SocCla)

p4 <- ggplot(SocCla, aes(x = SocCla, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") + 
  xlab("SocCla")+
  ylab("Count")+
  theme(axis.text.x = element_text(size=7,angle=0, color="skyblue4"))


ChildTrauma <- table[,7]
sum(is.na(ChildTrauma))
ChildTrauma <- na.omit(ChildTrauma)
ChildTrauma <- ChildTrauma %>% dplyr::count(ChildTrauma)

p5 <- ggplot(ChildTrauma, aes(x = ChildTrauma, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") + 
  xlab("ChildTrauma")+
  ylab("Count")+
  theme(axis.text.x = element_text(size=7,angle=0, color="skyblue4"))


Cannabis <- table[,8]
sum(is.na(Cannabis))
Cannabis <- na.omit(Cannabis)
Cannabis <- Cannabis %>% dplyr::count(Cannabis)

p6 <- ggplot(Cannabis, aes(x = Cannabis, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") + 
  xlab("Cannabis")+
  ylab("Count")+
  theme(axis.text.x = element_text(size=7,angle=0, color="skyblue4"))

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)


#Due to very few entries in some countries, lets classify them into the 'Other country' category (14)

#For the SocCla and EduLevel data ->
#levels(table$Country)
#plot(table$Country)
#table[,3] <- as.character(unlist(table[,3]))
#table$Country <- ifelse(table$Country == '13' | table$Country == '15' | table$Country == '3' | table$Country == '4' | 
#table$Country == '5' | table$Country == '79' | table$Country == '8', '14', table$Country)
#table[,3] <- as.factor(unlist(table[,3]))
#levels(table$Country)
#plot(table[,3])


#Correlation
pairs.panels(table[,-c(1,3:5,7:8)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             cex=5,
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)

#High correlation between PCA1 and ScoreInt and ScoreEA, thus, we should not consider PCA1


#Deleting Missing Values
table1 <- table[,-c(4:5)]
sum(is.na(table1))
table1 <- na.omit(table1)
#857 valid entries for the regression.

#Summary of the final data
summary(table1)
skim(table1)
head(table1)


#For the Income data ->
levels(table1$Country)
plot(table1$Country)
table1[,3] <- as.character(unlist(table1[,3]))
table1$Country <- ifelse(table1$Country == '11' | table1$Country == '12' | table1$Country == '13' | table1$Country == '3' | 
                           table1$Country == '4' | table1$Country == '5' | table1$Country == '6' | table1$Country == '79' | 
                           table1$Country == '8' | table1$Country == '9' | table1$Country == '15', '14', table1$Country)
table1[,3] <- as.factor(unlist(table1[,3]))
levels(table1$Country)
plot(table1[,3])


#Transformation
table1T=mutate(table1, Income=log(Income)^3)
hist(as.numeric(unlist(table1T[,4])))


########Standardization of Scores and PCAs 

table1.scaled <- scale(table1T[,7:19])
table2 <- cbind(table1.scaled, table1T[,-c(7:19)])
table2 <- table2[table2[,17]>0,]
#839 final entries


#Check that we get mean of 0 and sd of 1
colMeans(table1.scaled) 
apply(table1.scaled, 2, sd)


#################### ANALISIS OF EACH SCORE VARIABLE INDIVIDUALLY  ------------------------------

N <- length(table2[,17])
N
#839


#Null model
model0 <-  lm(Income ~ Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table2)
var0 <- summary(model0)$adj.r.squared
var0


#Explained variance for each of the Score variables

table_variances <- table2[,-c(1:10,14:19)]
vec <- colnames(table_variances)
R2 <- c(0)
pvalues <- c(0)
coefficients <- c(0)

for (i in 1:(length(vec))) {
  f <- as.formula(paste("Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10",vec[i],sep="+"))
  mod <- lm(f,table2)
  R2[i] <- summary(mod)$adj.r.squared
  pvalues[i] <- anova(model0,mod)[2,6]
  coefficients[i] <- summary(mod)$coeff[13]
}
variance_var <- data.frame(cbind(vec, R2, pvalues, coefficients))
variance_var[,c(2:4)] <- format(variance_var[,c(2:4)], decimal.mark = '.')
variance_var[,c(2:4)] <- as.numeric(unlist(variance_var[,c(2:4)]))
variance_var[2] <- variance_var[2]-var0 
variance_var[5] <- ifelse(variance_var[3]<0.1,"Significant","Not Significant") 
variance_var


#Significant: ScoreInt

plot(Effect(focal.predictors = "ScoreInt",lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table2)))




####################################### INTELIGENCE SCORES ---------------------------------------------------------------------


########### ANALISIS OF EACH ENVIRONMENTAL VARIABLE IN THE EFFECT OF SCOREINT ------------------------------


########### STRATIFICATION -------------------------------

#SEX
plot(table2[,14])
summary((table2[,14]))


#All Men=1
table3 <- table2[table2[,14]==1,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table3)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table3)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.008389977
anova(model0,model1)[2,6]
#P-value:0.009760275 Significant

coef(model1)[12]
#14.01431

plot(Effect(focal.predictors = "ScoreInt",model1))


#All Women = 2
table4 <- table2[table2[,14]==2,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table4)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table4)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.002233413
anova(model0,model1)[2,6]
#P-value:0.1635446 Not Significant

coef(model1)[12]
#7.332237

plot(Effect(focal.predictors = "ScoreInt",model1))



#CHILDTRAUMA
plot(table2[,18])
summary(table2[,18])

#All ChildTrauma=1
table5 <- table2[table2[,18]==1,]

model0 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table5)
model1 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table5)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.001728843
anova(model0,model1)[2,6]
#P-value:0.182464 Not significant

coef(model1)[13]
#8.968532

plot(Effect(focal.predictors = "ScoreInt",model1))


#All ChildTrauma = 0
table6 <- table2[table2[,18]==0,]

model0 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table6)
model1 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table6)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.01161277
anova(model0,model1)[2,6]
#P-value:0.004181987 Significant

coef(model1)[13]
#12.98921

plot(Effect(focal.predictors = "ScoreInt",model1))



#CANNABIS
plot(table2[,19])
summary(table2[,19])

#All Cannabis=1
table7 <- table2[table2[,19]==1,]

model0 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table7)
model1 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table7)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.004343122
anova(model0,model1)[2,6]
#P-value:0.1387727 Not significant

coef(model1)[13]
#11.54429

plot(Effect(focal.predictors = "ScoreInt",model1))


#All Cannabis = 0
table8 <- table2[table2[,19]==0,]

model0 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table8)
model1 <- lm(Income~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table8)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.005437038
anova(model0,model1)[2,6]
#P-value:0.02680654 Significant

coef(model1)[13]
#9.22918827

plot(Effect(focal.predictors = "ScoreInt",model1))




############## ANALYSIS OF SIGNIFICANT ENVIRONMENTAL VARIABLES IN RELATION TO GENDER AND COUNTRY IN THE EXPLAINED VARIANCE OF THE SCOREINT -----------------------

########### STRATIFICATION -------------------------------

#SEX

#All Men=1
table3 <- table2[table2[,14]==1,]
plot(table3[,19])
summary(table3[,19])

#If Cannabis=1
table9 <- table3[table3[,19]==1,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table9)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table9)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.02177849
anova(model0,model1)[2,6]
#P-value:0.01848771. Significant

coef(model1)[12]
#20.13763

plot(Effect(focal.predictors = "ScoreInt",model1))


#If Cannabis=0
table10 <- table3[table3[,19]==0,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table10)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table10)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:-0.0005781089
anova(model0,model1)[2,6]
#P-value:0.379758. Not Significant

coef(model1)[12]
#6.090072

plot(Effect(focal.predictors = "ScoreInt",model1))



#If ChildTrauma=1
plot(table3[,18])
summary(table3[,18])
table11 <- table3[table3[,18]==1,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table11)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table11)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.006012405
anova(model0,model1)[2,6]
#P-value:0.1060857 Not Significant

coef(model1)[12]
#14.63388

plot(Effect(focal.predictors = "ScoreInt",model1))


#If ChildTrauma=0
table12 <- table3[table3[,18]==0,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table12)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table12)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.004246341
anova(model0,model1)[2,6]
#P-value:0.1032806 Not Significant

coef(model1)[12]
#11.06117

plot(Effect(focal.predictors = "ScoreInt",model1))


#All Women=1
table4 <- table2[table2[,14]==2,]
plot(table4[,19])
summary(table4[,19])

#If Cannabis=1
table13 <- table4[table4[,19]==1,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table13)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table13)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:-0.001363904
anova(model0,model1)[2,6]
#P-value:0.3415064 Not Significant

coef(model1)[12]
#-17.94864

plot(Effect(focal.predictors = "ScoreInt",model1))


#If Cannabis=0
table14 <- table4[table4[,19]==0,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table14)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table14)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.005459589
anova(model0,model1)[2,6]
#P-value:0.09796774 Significant

coef(model1)[12]
#8.80458

plot(Effect(focal.predictors = "ScoreInt",model1))



#If ChildTrauma=1
plot(table4[,18])
summary(table4[,18])
table15 <- table4[table4[,18]==1,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table15)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table15)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:-0.004944117
anova(model0,model1)[2,6]
#P-value:0.6853207 Not Significant

coef(model1)[12]
#-4.370429

plot(Effect(focal.predictors = "ScoreInt",model1))


#If ChildTrauma=0
table16 <- table4[table4[,18]==0,]

model0 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table16)
model1 <- lm(Income~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table16)
summary(model0)
summary(model1)

var0 <- summary(model0)$adj.r.squared
var1 <- summary(model1)$adj.r.squared
var1-var0
#R2:0.01218105
anova(model0,model1)[2,6]
#P-value:0.04976168 Significant

coef(model1)[12]
#12.55072

plot(Effect(focal.predictors = "ScoreInt",model1))



#COUNTRY

#If Men and Cannabis=1
table9 <- table3[table3[,19]==1,]
plot(table9[,16])
#Very few observations to study the differences in both countries




#############################  STEPWISE BEST MODEL SELECTION -----------------------------------------------------------

#Models with only genetic variables

# selection with AIC criterion:
set.seed(3890)
fit.all <- lm(Income ~ Age+Sex+ScoreInt+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data = table2)
model0 = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
model0 <- lm(Income ~ Age + Sex + ScoreInt + PCA2 + PCA3 + PCA4 + 
               PCA8 + PCA10 + PCA2:PCA10 + PCA4:PCA10 + Age:PCA2 + Age:Sex, data = table2)
extractAIC(model0)[2] 
#7624.151
summary(model0)$adj.r.squared
#0.3346525

coef(model0)[4]
#10.46875


#Models with environmental factors 

# selection with AIC criterion:
set.seed(3890)
fit.all <- lm(Income ~ Age+Sex+ScoreInt+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+Country+ChildTrauma+Cannabis, data = table2)
best_model = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
best_model <- lm(Income ~ Age + Sex + ScoreInt + PCA8 + PCA10 + Country + 
                   Cannabis + PCA10:Country + Age:Cannabis + PCA8:Cannabis + 
                   Age:Sex, data = table2)
extractAIC(best_model)[2] 
#7569.096
summary(best_model)$adj.r.squared
#0.3776413

coef(best_model)[4]
#10.93616



#How much the score impact in the incomes:
scoremean <- mean(table1$ScoreInt)
scoresd <- sd(table1$ScoreInt)
sc = destandardized(1,scoremean,scoresd)
#1score=0.0002083181

#Des-transformation
#1euro=log(2)^3=0.3330247

#Equivalence:
#1 unit more of score -> 10.94 euros
#0.0002083181 -> 32.85 euros
#1 -> 157,697.7 euros
#0.00001 -> 1.58 euros



##################### EVALUATION OF THE BEST MODEL -------------------------------------------------------

#Plot of the model
plot(best_model)
# 1: Linearity: good because the line is horizontal
# 2: Normal distribution: good
# 3: Homoscedasticity: more or less constant variance of errors



#########################  PLOTING THE EFFECTS OF THE VARIABLES IN THE BEST MODEL ------------------------------------------------------

plot(Effect(focal.predictors = "Age",best_model))
#Older people have more probability to have higher Incomes

plot(Effect(focal.predictors = "Sex",best_model))
#The women have much lower probability to have higher Incomes

plot(Effect(focal.predictors = "Country",best_model))
#The higher the Score, the higher the Incomes

plot(Effect(focal.predictors = "ScoreInt",best_model))
#The higher the Score, the higher the Incomes

plot(Effect(focal.predictors = "Cannabis",best_model))
#The consumption of Cannabis involves to have less probability to gain higher Incomes



#########
## REGRESSION GWAS SOCIAL CLASS
########


rm(list=ls())



###########################################  PREPROCESSING -----------------------------------------------------------

table <- table[,-c(1:2)]
table <- table[table[,1]==3,]
table <- table[,-c(1,5,7)]


names(table)[4] <- c("SocCla")
names(table)[5] <- c("ChildTrauma")
names(table)[6] <- c("Cannabis")
names(table)[17] <- c("ScoreInt")
names(table)[18] <- c("ScoreCog")
names(table)[19] <- c("ScoreEA")

table[table=="NA"] <- NA

table[,1] <- as.factor(unlist((table[,1])))
table[,2] <- as.numeric(unlist(table[,2]))
table[,3] <- as.factor(unlist(table[,3]))
table[,4] <- as.factor(unlist(table[,4]))
table[,5] <- as.factor(unlist(table[,5]))
table[,6] <- as.factor(unlist(table[,6]))


for (i in 7:19) {
  table[,i] <- as.numeric(unlist(table[,i]))
}

#Deleting Missing Values
sum(is.na(table))
table1 <- na.omit(table)
#2,309 valid entries for the regression.


#Grouping the Phenotypes in only 3 groups
table1[,4] <- as.character(unlist(table1[,4]))

table1$SocCla <- ifelse(table1$SocCla == '1' | table1$SocCla == '2', '1', table1$SocCla)
table1$SocCla <- ifelse(table1$SocCla == '3' | table1$SocCla == '4', '2', table1$SocCla)
table1$SocCla <- ifelse(table1$SocCla == '5' | table1$SocCla == '6', '3', table1$SocCla)

table1[,4] <- as.factor(unlist(table1[,4]))
summary(table1$SocCla)

#Grouping the Countries with less observations
#Due to very low entries in some countries, lets classify them into the 'Other country' category (14)

levels(table1$Country)
plot(table1$Country)
table1[,3] <- as.character(unlist(table1[,3]))
table1$Country <- ifelse(table1$Country == '13' | table1$Country == '15' | table1$Country == '3' | table1$Country == '4' | 
                           table1$Country == '5' | table1$Country == '79' | table1$Country == '8', '14', table1$Country)
table1[,3] <- as.factor(unlist(table1[,3]))
levels(table1$Country)
plot(table1[,3])


########Standardization of Scores and PCAs 

table1.scaled <- scale(table1[,7:19])
table2 <- cbind(table1.scaled, table1[,-c(7:19)])

#Check that we get mean of 0 and sd of 1
colMeans(table1.scaled) 
apply(table1.scaled, 2, sd)


#################### ANALISIS OF EACH SCORE VARIABLE INDIVIDUALLY  ------------------------------

#Checking the proportional odds assumptions
sf <- function(y) {
  c('Y>=1' = qlogis(mean(y >= 1)),
    'Y>=2' = qlogis(mean(y >= 2)),
    'Y>=3' = qlogis(mean(y >= 3)))
}

s <- with(table2, summary(as.numeric(SocCla) ~ Age+Sex+Country+ScoreInt+ScoreCog+ScoreEA+ChildTrauma+Cannabis+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, fun=sf))
s[, 4] <- s[, 4] - s[, 3]
s[, 3] <- s[, 3] - s[, 3]
s

#The proportional odds don't meet completely between the categories of Age and PCA2 variables. That could be a problem in the efficiency of the
#models

N <- length(table2[,17])
N
#2.309


#Null model
model_null <- polr(SocCla~1,data=table2, Hess=T, method='logistic')
model0 <-  polr(SocCla ~ Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table2, Hess=T, method='logistic')
var0 <- 1-(model0$deviance/model_null$deviance)
var0


#Explained variance for each of the Score variables

table_variances <- table2[,c(11:13)]
vec <- colnames(table_variances)
R2 <- c(0)
coefficients <- c(0)
anovatest <- c(0)

for (i in 1:(length(vec))) {
  f <- as.formula(paste("SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10",vec[i],sep="+"))
  mod <- polr(f,table2, Hess=T, method='logistic')
  R2[i] <- 1-(mod$deviance/model0$deviance)
  coefficients[i] <- summary(mod)$coeff[12]
  anovatest[i] <- anova(model0,mod)[2,7]
}
variance_var <- data.frame(cbind(vec, R2, coefficients, anovatest))
variance_var[,c(2:4)] <- format(variance_var[,c(2:4)], decimal.mark = '.')
variance_var[,c(2:4)] <- as.numeric(unlist(variance_var[,c(2:4)]))
variance_var[5] <- ifelse(variance_var[4]<0.1,"Significant","Not Significant") 
variance_var



#Significant: ScoreInt and ScoreEA

plot(Effect(focal.predictors = "ScoreInt",polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table2, Hess=T, method='logistic')))
plot(Effect(focal.predictors = "ScoreEA",polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table2, Hess=T, method='logistic')))


####################################### INTELIGENCE SCORES ---------------------------------------------------------------------


########### ANALISIS OF EACH ENVIRONMENTAL VARIABLE IN THE EFFECT OF SCOREINT ------------------------------


########### STRATIFICATION -------------------------------

#SEX
plot(table2[,14])
summary((table2[,14]))

#All Men=1
table3 <- table2[table2[,14]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table3,  Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table3,  Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0008788934
anova(model0,model1)[2,7]
#P-value:0.1820796 Not significant

coef(model1)[11]
#0.08169811

plot(Effect(focal.predictors = "ScoreInt",model1))


#All Women = 2
table4 <- table2[table2[,14]==2,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table4, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table4, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002340389
anova(model0,model1)[2,7]
#P-value:0.01652412 Significant

coef(model1)[11]
#0.1293975

plot(Effect(focal.predictors = "ScoreInt",model1))



#CHILDTRAUMA
plot(table3[,18])
summary((table3[,18]))

#All ChildTrauma=1
table5 <- table2[table2[,18]==1,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table5, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table5, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.000106304
anova(model0,model1)[2,7]
#P-value:0.629511 Not significant

coef(model1)[12]
#0.02723482

plot(Effect(focal.predictors = "ScoreInt",model1))


#All ChildTrauma = 0
table6 <- table2[table2[,18]==0,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table6, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table6, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004571533
anova(model0,model1)[2,7]
#P-value:0.001000891. Significant

coef(model1)[12]
#0.1896252

plot(Effect(focal.predictors = "ScoreInt",model1))



#CANNABIS
plot(table2[,19])
summary((table2[,19]))

#All Cannabis=1
table7 <- table2[table2[,19]==1,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table7, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table7, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0004467137
anova(model0,model1)[2,7]
#P-value:0.390768 Not significant

coef(model1)[12]
#0.05453648

plot(Effect(focal.predictors = "ScoreInt",model1))


#All Cannabis = 0
table8 <- table2[table2[,19]==0,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table8,  Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table8, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002494812
anova(model0,model1)[2,7]
#P-value:0.007714229 Significant

coef(model1)[12]
#0.1457672

plot(Effect(focal.predictors = "ScoreInt",model1))




############## ANALYSIS OF SIGNIFICANT ENVIRONMENTAL VARIABLES IN RELATION TO GENDER AND COUNTRY IN THE EXPLAINED VARIANCE OF THE SCOREINT -----------------------

########### STRATIFICATION -------------------------------

#SEX

#All Men=1
table3 <- table2[table2[,14]==1,]
plot(table3[,19])
summary((table3[,19]))

#If Cannabis=1
table9 <- table3[table3[,19]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table9, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table9, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.000694219
anova(model0,model1)[2,7]
#P-value:0.410886 Not Significant

coef(model1)[11]
#0.07220459

plot(Effect(focal.predictors = "ScoreInt",model1))


#If Cannabis=0
table10 <- table3[table3[,19]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table10,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table10,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0003531614
anova(model0,model1)[2,7]
#P-value:0.5465277 Not Significant

coef(model1)[11]
#0.05398283

plot(Effect(focal.predictors = "ScoreInt",model1))



#If ChildTrauma=1
plot(table3[,18])
summary((table3[,18]))
table11 <- table3[table3[,18]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table11,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table11,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0006217189
anova(model0,model1)[2,7]
#P-value:0.4418533 Not Significant

coef(model1)[11]
#0.06640069

plot(Effect(focal.predictors = "ScoreInt",model1))


#If ChildTrauma=0
table12 <- table3[table3[,18]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table12,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table12,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001590386
anova(model0,model1)[2,7]
#P-value:0.1955179 Significant

coef(model1)[11]
#0.1151865

plot(Effect(focal.predictors = "ScoreInt",model1))



#All Women=2
table4 <- table2[table2[,14]==2,]
plot(table4[,19])
summary((table4[,19]))

#If Cannabis=1
table13 <- table4[table4[,19]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table13, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table13, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0002715093
anova(model0,model1)[2,7]
#P-value:0.6729238 Not Significant

coef(model1)[11]
#0.03973708

plot(Effect(focal.predictors = "ScoreInt",model1))


#If Cannabis=0
table14 <- table4[table4[,19]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table14,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table14,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004605897
anova(model0,model1)[2,7]
#P-value:0.004498022 Significant

coef(model1)[11]
#0.2019167

plot(Effect(focal.predictors = "ScoreInt",model1))



#If ChildTrauma=1
plot(table4[,18])
summary((table4[,18]))
table15 <- table4[table4[,18]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table15,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table15,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:2.128571e-05
anova(model0,model1)[2,7]
#P-value:0.8737156 Not Significant

coef(model1)[11]
#0.01177709

plot(Effect(focal.predictors = "ScoreInt",model1))


#If ChildTrauma=0
table16 <- table4[table4[,18]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table16,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table16,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009823552
anova(model0,model1)[2,7]
#P-value:0.0004475407 Significant

coef(model1)[11]
#0.2856033

plot(Effect(focal.predictors = "ScoreInt",model1))



#COUNTRY


#If Men and Cannabis=0
table10 <- table3[table3[,19]==0,]
plot(table10[,16]) #Lets analyze only Turkey, Brazil and Spain
summary(table10[,16])

#Turkey
table17 <- table10[table10[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table17,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table17,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.006202
anova(model0,model1)[2,7]
#P-value:0.1139535 Not Significant

coef(model1)[11]
#-0.7744936

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table18 <- table10[table10[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table18,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table18,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002782281
anova(model0,model1)[2,7]
#P-value:0.5302118 Not Significant

coef(model1)[11]
#-0.2037396

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table19 <- table10[table10[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table19,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table19,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01815179
anova(model0,model1)[2,7]
#P-value:0.144781 Not Significant

coef(model1)[11]
#-1.585723

plot(Effect(focal.predictors = "ScoreInt",model1))


#Others
table29 <- table10[table10[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table29,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table29,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009042275
anova(model0,model1)[2,7]
#P-value:0.2847705 Not Significant

coef(model1)[11]
#0.1560957

plot(Effect(focal.predictors = "ScoreInt",model1))



#If Men and ChildTrauma=0
table12 <- table3[table3[,18]==0,]
plot(table12[,16])
summary(table12[,16])


#Turkey
table20 <- table12[table12[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table20,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table20,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.005424876
anova(model0,model1)[2,7]
#P-value:0.3510999 Not Significant

coef(model1)[11]
#-0.7410346

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table21 <- table12[table12[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table21,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table21,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0001075477
anova(model0,model1)[2,7]
#P-value:0.9073118 Significant

coef(model1)[11]
#0.044034

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table22 <- table12[table12[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table22,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table22,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:9.371988e-05
anova(model0,model1)[2,7]
#P-value:0.8521723 Not Significant

coef(model1)[11]
#0.08712374

plot(Effect(focal.predictors = "ScoreInt",model1))


#Others
table30 <- table12[table12[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table30,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table30,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002453971
anova(model0,model1)[2,7]
#P-value:0.5816736 Not Significant

coef(model1)[11]
#0.09309562

plot(Effect(focal.predictors = "ScoreInt",model1))



#If Women and Cannabis=0
table14 <- table4[table4[,19]==0,]
plot(table14[,16])
summary(table14[,16])


#Turkey
table23 <- table14[table14[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table23,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table23,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.007090144
anova(model0,model1)[2,7]
#P-value:0.02537401 Not Significant

coef(model1)[11]
#0.8066957

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table24 <- table14[table14[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table24,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table24,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0269518
anova(model0,model1)[2,7]
#P-value:0.011678 Significant

coef(model1)[11]
#0.6882456

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table25 <- table14[table14[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table25,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table25,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.03799573
anova(model0,model1)[2,7]
#P-value:0.01166374 Significant

coef(model1)[11]
#1.781989

plot(Effect(focal.predictors = "ScoreInt",model1))


#Others
table31 <- table14[table14[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table31,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table31,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0002727653
anova(model0,model1)[2,7]
#P-value:0.8087427 Significant

coef(model1)[11]
#0.03226379

plot(Effect(focal.predictors = "ScoreInt",model1))



#If Women and ChildTrauma=0
table16 <- table4[table4[,18]==0,]
plot(table16[,16])
summary(table16[,16])


#Turkey
table26 <- table16[table16[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table26,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table26,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01345715
anova(model0,model1)[2,7]
#P-value:0.0542441 Significant

coef(model1)[11]
#1.194584

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table27 <- table16[table16[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table27,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table27,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.07918945
anova(model0,model1)[2,7]
#P-value:0.0004363481 Significant

coef(model1)[11]
#1.115877

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table28 <- table16[table16[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table28,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table28,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:9.665973e-06
anova(model0,model1)[2,7]
#P-value:0.9610299 Not Significant

coef(model1)[11]
#0.02847335

plot(Effect(focal.predictors = "ScoreInt",model1))


#Others
table32 <- table16[table16[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table32,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table32,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001737059
anova(model0,model1)[2,7]
#P-value:0.5727678 Not Significant

coef(model1)[11]
#0.1034365

plot(Effect(focal.predictors = "ScoreInt",model1))




#############################  STEPWISE BEST MODEL SELECTION -----------------------------------------------------------

#Models with only genetic variables

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(SocCla ~ Age+Sex+ScoreInt+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data = table2, Hess=T, method='logistic')
model0 = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
model0 <- polr(SocCla ~ Age + Sex + ScoreInt + PCA2 + PCA3 + 
                 PCA4 + PCA8 + PCA9 + Sex:PCA2 + Age:Sex + Age:PCA4 + Age:PCA2 + 
                 Sex:PCA3 + PCA8:PCA9 + ScoreInt:PCA4 + Sex:PCA9, data = table2, Hess=T, method='logistic')
extractAIC(model0)[2] 
#4471.969
var1 <- 1-(model0$deviance/model_null$deviance)
var1
#R2:0.06699514

coef(model0)[3]
#0.1894415


#Models with environmental factors 

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(SocCla ~ Age+Sex+ScoreInt+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+Country+ChildTrauma+Cannabis, data = table2, Hess=T, method='logistic')
best_model = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
best_model <- polr(SocCla ~ Age + Sex + ScoreInt + PCA5 + PCA8 + 
                     PCA9 + Country + ChildTrauma + Sex:Country + Age:Country + 
                     Age:Sex + PCA8:PCA9 + Age:PCA5 + ScoreInt:ChildTrauma + Age:ScoreInt, data = table2, Hess=T, method='logistic')
extractAIC(best_model)[2] 
#4407.671
var1 <- 1-(best_model$deviance/model_null$deviance)
var1
#R2:0.08640796

coef(best_model)[3]
#0.5020675



####################################### EDUCATIONAL ATTAINMENT SCORES ---------------------------------------------------------------------


########### ANALISIS OF EACH ENVIRONMENTAL VARIABLE IN THE EFFECT OF SCOREINT ------------------------------


########### STRATIFICATION -------------------------------

#SEX
plot(table2[,14])
summary(table2[,14])

#All Men=1
table3 <- table2[table2[,14]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table3,  Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table3,  Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001274298
anova(model0,model1)[2,7]
#P-value:0.1081105 Not significant

coef(model1)[11]
#0.09672878

plot(Effect(focal.predictors = "ScoreEA",model1))


#All Women = 2
table4 <- table2[table2[,14]==2,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table4, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table4, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.003346177
anova(model0,model1)[2,7]
#P-value:0.004152997 Significant

coef(model1)[11]
#0.156337

plot(Effect(focal.predictors = "ScoreEA",model1))



#CHILDTRAUMA
plot(table3[,18])

#All ChildTrauma=1
table5 <- table2[table2[,18]==1,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table5, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table5, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0002359879
anova(model0,model1)[2,7]
#P-value:0.4722814 Not significant

coef(model1)[12]
#0.04013365

plot(Effect(focal.predictors = "ScoreEA",model1))


#All ChildTrauma = 0
table6 <- table2[table2[,18]==0,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table6, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table6, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.005720564
anova(model0,model1)[2,7]
#P-value:0.0002326686 Significant

coef(model1)[12]
#0.2139361

plot(Effect(focal.predictors = "ScoreEA",model1))



#CANNABIS
plot(table2[,19])

#All Cannabis=1
table7 <- table2[table2[,19]==1,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table7, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table7, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:9.855517e-05
anova(model0,model1)[2,7]
#P-value:0.6868653 Not significant

coef(model1)[12]
#0.02585325

plot(Effect(focal.predictors = "ScoreEA",model1))


#All Cannabis = 0
table8 <- table2[table2[,19]==0,]

model0 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table8,  Hess=T, method='logistic')
model1 <- polr(SocCla~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table8, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004248555
anova(model0,model1)[2,7]
#P-value:0.0005072895 Significant

coef(model1)[12]
#0.1886352

plot(Effect(focal.predictors = "ScoreEA",model1))




############## ANALYSIS OF SIGNIFICANT ENVIRONMENTAL VARIABLES IN RELATION TO GENDER AND COUNTRY IN THE EXPLAINED VARIANCE OF THE ScoreEA -----------------------

########### STRATIFICATION -------------------------------

#SEX

#All Men=1
table3 <- table2[table2[,14]==1,]
plot(table3[,19])

#If Cannabis=1
table9 <- table3[table3[,19]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table9, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table9, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:7.75904e-05
anova(model0,model1)[2,7]
#P-value:0.7833785 Not Significant

coef(model1)[11]
#0.02454292

plot(Effect(focal.predictors = "ScoreEA",model1))


#If Cannabis=0
table10 <- table3[table3[,19]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table10,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table10,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002232239
anova(model0,model1)[2,7]
#P-value:0.1295359 Not Significant

coef(model1)[11]
#0.1276956

plot(Effect(focal.predictors = "ScoreEA",model1))



#If ChildTrauma=1
table11 <- table3[table3[,18]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table11,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table11,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0003020007
anova(model0,model1)[2,7]
#P-value:0.5919529 Not Significant

coef(model1)[11]
#0.04535266

plot(Effect(focal.predictors = "ScoreEA",model1))


#If ChildTrauma=0
table12 <- table3[table3[,18]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table12,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table12,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.003563223
anova(model0,model1)[2,7]
#P-value:0.05268093 Significant

coef(model1)[11]
#0.1692685

plot(Effect(focal.predictors = "ScoreEA",model1))



#All Women=2
table4 <- table2[table2[,14]==2,]
plot(table4[,19])

#If Cannabis=1
table13 <- table4[table4[,19]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table13, Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table13, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0003994911
anova(model0,model1)[2,7]
#P-value:0.6086127 Not Significant

coef(model1)[11]
#0.04812528

plot(Effect(focal.predictors = "ScoreEA",model1))


#If Cannabis=0
table14 <- table4[table4[,19]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table14,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table14,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.006191401
anova(model0,model1)[2,7]
#P-value:0.0009883541 Significant

coef(model1)[11]
#0.2373782

plot(Effect(focal.predictors = "ScoreEA",model1))



#If ChildTrauma=1
summary(table4[,18])
table15 <- table4[table4[,18]==1,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table15,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table15,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0006028783
anova(model0,model1)[2,7]
#P-value:0.3976227 Not Significant

coef(model1)[11]
#0.06314846

plot(Effect(focal.predictors = "ScoreEA",model1))


#If ChildTrauma=0
table16 <- table4[table4[,18]==0,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table16,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table16,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009837817
anova(model0,model1)[2,7]
#P-value:0.0004432711 Significant

coef(model1)[11]
#0.2856033

plot(Effect(focal.predictors = "ScoreEA",model1))



#COUNTRY


#If Men and Cannabis=0
table10 <- table3[table3[,19]==0,]
plot(table10[,16]) #Lets analyze only Turkey, Brazil and Spain

#Turkey
table17 <- table10[table10[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table17,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table17,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0001002269
anova(model0,model1)[2,7]
#P-value:0.8407451 Not Significant

coef(model1)[11]
#-0.07358145

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table18 <- table10[table10[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table18,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table18,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0001213244
anova(model0,model1)[2,7]
#P-value:0.8957174 Not Significant

coef(model1)[11]
#0.03973474

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table19 <- table10[table10[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table19,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table19,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01812602
anova(model0,model1)[2,7]
#P-value:0.1450665 Not Significant

coef(model1)[11]
#0.9036463

plot(Effect(focal.predictors = "ScoreEA",model1))


#Others
table29 <- table10[table10[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table29,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table29,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.008863697
anova(model0,model1)[2,7]
#P-value:0.2895775 Not Significant

coef(model1)[11]
#0.1560899

plot(Effect(focal.predictors = "ScoreEA",model1))



#If Men and ChildTrauma=0
table12 <- table3[table3[,18]==0,]
plot(table12[,16])


#Turkey
table20 <- table12[table12[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table20,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table20,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.000526953
anova(model0,model1)[2,7]
#P-value:0.7713446 Not Significant

coef(model1)[11]
#-0.1921871

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table21 <- table12[table12[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table21,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table21,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.006643406
anova(model0,model1)[2,7]
#P-value:0.3601488 Significant

coef(model1)[11]
#0.3118827

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table22 <- table12[table12[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table22,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table22,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.00135199
anova(model0,model1)[2,7]
#P-value:0.4790858 Not Significant

coef(model1)[11]
#0.2092054

plot(Effect(focal.predictors = "ScoreEA",model1))


#Others
table30 <- table12[table12[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table30,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table30,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0008461644
anova(model0,model1)[2,7]
#P-value:0.746303 Not Significant

coef(model1)[11]
#0.06153446

plot(Effect(focal.predictors = "ScoreEA",model1))



#If Women and Cannabis=0
table14 <- table4[table4[,19]==0,]
plot(table14[,16])


#Turkey
table23 <- table14[table14[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table23,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table23,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.007107832
anova(model0,model1)[2,7]
#P-value:0.02519188 Not Significant

coef(model1)[11]
#0.6524434

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table24 <- table14[table14[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table24,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table24,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.02445884
anova(model0,model1)[2,7]
#P-value:0.01629357 Significant

coef(model1)[11]
#0.6146912

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table25 <- table14[table14[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table25,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table25,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.04890483
anova(model0,model1)[2,7]
#P-value:0.004217557 Significant

coef(model1)[11]
#1.541436

plot(Effect(focal.predictors = "ScoreEA",model1))


#Others
table31 <- table14[table14[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table31,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table31,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:9.792497e-05
anova(model0,model1)[2,7]
#P-value:0.8846881 Significant

coef(model1)[11]
#-0.01943644

plot(Effect(focal.predictors = "ScoreEA",model1))



#If Women and ChildTrauma=0
table16 <- table4[table4[,18]==0,]
plot(table16[,16])


#Turkey
table26 <- table16[table16[,16]==10,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table26,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table26,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0004098672
anova(model0,model1)[2,7]
#P-value:0.7369234 Not Significant

coef(model1)[11]
#0.1620919

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table27 <- table16[table16[,16]==12,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table27,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table27,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.04175249
anova(model0,model1)[2,7]
#P-value:0.01065527 Significant

coef(model1)[11]
#0.7896875

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table28 <- table16[table16[,16]==7,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table28,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table28,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.005220567
anova(model0,model1)[2,7]
#P-value:0.2561514 Not Significant

coef(model1)[11]
#0.5227248

plot(Effect(focal.predictors = "ScoreEA",model1))


#Others
table32 <- table16[table16[,16]==14,]

model0 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table32,Hess=T, method='logistic')
model1 <- polr(SocCla~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table32,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0005103323
anova(model0,model1)[2,7]
#P-value:0.7598398 Not Significant

coef(model1)[11]
#0.05589382

plot(Effect(focal.predictors = "ScoreEA",model1))





#############################  STEPWISE BEST MODEL SELECTION -----------------------------------------------------------

#Models with only genetic variables

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(SocCla ~ Age+Sex+ScoreEA+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data = table2, Hess=T, method='logistic')
model0 = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
model0 <- polr(SocCla ~ Age + Sex + ScoreEA + PCA2 + PCA3 + PCA4 + 
                 PCA5 + PCA6 + PCA8 + PCA9 + PCA10 + Sex:PCA2 + Age:Sex + 
                 Age:PCA4 + Age:PCA2 + Sex:PCA3 + PCA8:PCA9 + ScoreEA:PCA6 + 
                 Sex:PCA5 + PCA4:PCA10 + Sex:PCA9, data = table2, Hess=T, method='logistic')
extractAIC(model0)[2] 
#4461.848
var0 <- 1-(model0$deviance/model_null$deviance)
var0
#0.07122722

coef(model0)[3]
#0.1318695


#Models with environmental factors 

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(SocCla ~ Age+Sex+ScoreEA+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+Country+ChildTrauma+Cannabis, data = table2, Hess=T, method='logistic')
best_model = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
best_model <- polr(SocCla ~ Age + Sex + ScoreEA + PCA5 + PCA6 + PCA8 + 
                     PCA9 + Country + ChildTrauma + Sex:Country + Age:Country + 
                     Age:Sex + PCA8:PCA9 + Age:PCA5 + ScoreEA:PCA6 + ScoreEA:ChildTrauma, data = table2, Hess=T, method='logistic')
extractAIC(best_model)[2] 
#4401.623
var0 <- 1-(best_model$deviance/model_null$deviance)
var0
#0.08810075

coef(best_model)[3]
#0.2534786

#Best model is with ScoreEA



#########################  PLOTING THE EFFECTS OF THE VARIABLES IN THE BEST MODEL ------------------------------------------------------

plot(Effect(focal.predictors = "Age",best_model))
#Older people have more probability to have higher SocClas

plot(Effect(focal.predictors = "Sex",best_model))
#The women have much lower probability to have higher SocClas

plot(Effect(focal.predictors = "Country",best_model))
#The higher the Score, the higher the SocClas

plot(Effect(focal.predictors = "ScoreEA",best_model))
#The higher the Score, the higher the SocClas

plot(Effect(focal.predictors = "ChildTrauma",best_model))
#The consumption of ChildTrauma involves to have less probability to gain higher SocClas

plot(Effect(focal.predictors = "Cannabis",best_model))
#Not significant



################# PREDICTION POWER OF THE BEST MODEL -------------------------------------------

#SENSIBILITY TEST
ROC(form=SocCla ~ Age + Sex + ScoreEA + PCA5 + PCA6 + PCA8 + 
      PCA9 + Country + ChildTrauma + Sex:Country + Age:Country + 
      Age:Sex + PCA8:PCA9 + Age:PCA5 + ScoreEA:PCA6 + ScoreEA:ChildTrauma, data=table2, plot="ROC", lwd=3, cex=1.5, AUC=TRUE)
#AUC=0.790
#Acceptable AUC for classification with 3 levels



#CONFUSION MATRIX

#Prepare Training and Test Data
set.seed(3890)
trainingRows <- sample(1:nrow(table2), 0.7 * nrow(table2))
trainingData <- table2[trainingRows, ]
testData <- table2[-trainingRows, ]

finalmod2 <- polr(SocCla ~ Age + Sex + ScoreEA + PCA5 + PCA6 + PCA8 + 
                    PCA9 + Country + ChildTrauma + Sex:Country + Age:Country + 
                    Age:Sex + PCA8:PCA9 + Age:PCA5 + ScoreEA:PCA6 + ScoreEA:ChildTrauma, data=table2, Hess=TRUE, method="logistic")
summary(finalmod2)

predictedPhenotype <- predict(finalmod2, testData)  # predict the classes directly
head(predictedPhenotype)

predictedScores <- predict(finalmod2, testData, type="p")  # predict the probabilites
head(predictedScores)

ConfMat <- table(testData$SocCla, predictedPhenotype)  # confusion matrix
ConfMat

mean(as.character(testData$SocCla) != as.character(predictedPhenotype))  # missclassification error
#A 42.2% of error.


n = length(testData$SocCla)
accuracy <- sum(diag(ConfMat)) / n
accuracy
#A 57.7% of accuracy.










#########
## REGRESSION GWAS EDUCATION LEVEL
########

rm(list=ls())


###########################################  PREPROCESSING -----------------------------------------------------------

table <- table[,-c(1:2)]
table <- table[table[,1]==3,]
table <- table[,-c(1,6,7)]

names(table)[4] <- c("EduLevel")
names(table)[5] <- c("ChildTrauma")
names(table)[6] <- c("Cannabis")
names(table)[17] <- c("ScoreInt")
names(table)[18] <- c("ScoreCog")
names(table)[19] <- c("ScoreEA")

table[table=="NA"] <- NA

table[,1] <- as.factor(unlist((table[,1])))
table[,2] <- as.numeric(unlist(table[,2]))
table[,3] <- as.factor(unlist(table[,3]))
table[,4] <- as.factor(unlist(table[,4]))
table[,5] <- as.factor(unlist(table[,5]))
table[,6] <- as.factor(unlist(table[,6]))


for (i in 7:19) {
  table[,i] <- as.numeric(unlist(table[,i]))
}

#Deleting Missing Values
sum(is.na(table))
table1 <- na.omit(table)
#2,301 valid entries for the regression.

#Grouping the Phenotypes in only 3 groups
table1[,4] <- as.numeric(unlist(table1[,4]))

table1$EduLevel <- ifelse(table1$EduLevel == 1 | table1$EduLevel == 2, '1', table1$EduLevel)
table1$EduLevel <- ifelse(table1$EduLevel == 3 | table1$EduLevel == 4, '2', table1$EduLevel)
table1$EduLevel <- ifelse(table1$EduLevel == 5 | table1$EduLevel == 6, '3', table1$EduLevel)

table1[,4] <- as.factor(unlist(table1[,4]))


#Grouping the Countries with less observations
#Due to very low entries in some countries, lets classify them into the 'Other country' category (14)

levels(table1$Country)
plot(table1$Country)
table1[,3] <- as.character(unlist(table1[,3]))
table1$Country <- ifelse(table1$Country == '13' | table1$Country == '15' | table1$Country == '3' | table1$Country == '4' | 
                           table1$Country == '5' | table1$Country == '79' | table1$Country == '8', '14', table1$Country)
table1[,3] <- as.factor(unlist(table1[,3]))
levels(table1$Country)
plot(table1[,3])

########Standardization of Scores and PCAs 

table1.scaled <- scale(table1[,7:19])
table2 <- cbind(table1.scaled, table1[,-c(7:19)])

#Check that we get mean of 0 and sd of 1
colMeans(table1.scaled) 
apply(table1.scaled, 2, sd)




############ INTELIGENCE SCORES #########################################

#Checking the proportional odds assumptions
sf <- function(y) {
  c('Y>=1' = qlogis(mean(y >= 1)),
    'Y>=2' = qlogis(mean(y >= 2)),
    'Y>=3' = qlogis(mean(y >= 3)))
}

s <- with(table2, summary(as.numeric(EduLevel) ~ Age+Sex+Country+ScoreInt+ScoreCog+ScoreEA+ChildTrauma+Cannabis+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, fun=sf))
s[, 4] <- s[, 4] - s[, 3]
s[, 3] <- s[, 3] - s[, 3]
s

#The proportional odds don't meet completely between the categories of Age and PCA2 variables. That could be a problem in the efficiency of the
#models

N <- length(table2[,17])
N
#2.301

#Null model
model_null <- polr(EduLevel~1,data=table2, Hess=T, method='logistic')
model0 <-  polr(EduLevel ~ Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table2, Hess=T, method='logistic')
var0 <- 1-(model0$deviance/model_null$deviance)
var0


#Explained variance for each of the Score variables

table_variances <- table2[,c(11:13)]
vec <- colnames(table_variances)
R2 <- c(0)
coefficients <- c(0)
anovatest <- c(0)

for (i in 1:(length(vec))) {
  f <- as.formula(paste("EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10",vec[i],sep="+"))
  mod <- polr(f,table2, Hess=T, method='logistic')
  R2[i] <- 1-(mod$deviance/model0$deviance)
  coefficients[i] <- summary(mod)$coeff[12]
  anovatest[i] <- anova(model0,mod)[2,7]
}
variance_var <- data.frame(cbind(vec, R2, coefficients, anovatest))
variance_var[,c(2:4)] <- format(variance_var[,c(2:4)], decimal.mark = '.')
variance_var[,c(2:4)] <- as.numeric(unlist(variance_var[,c(2:4)]))
variance_var[5] <- ifelse(variance_var[4]<0.1,"Significant","Not Significant") 
variance_var



#Significant: ScoreInt and ScoreEA

plot(Effect(focal.predictors = "ScoreInt",polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table2, Hess=T, method='logistic')))
plot(Effect(focal.predictors = "ScoreCog",polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table2, Hess=T, method='logistic')))
plot(Effect(focal.predictors = "ScoreEA",polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table2, Hess=T, method='logistic')))


####################################### INTELIGENCE SCORES ---------------------------------------------------------------------


########### ANALISIS OF EACH ENVIRONMENTAL VARIABLE IN THE EFFECT OF SCOREINT ------------------------------


########### STRATIFICATION -------------------------------

#SEX
plot(table2[,14])
summary(table2[,14])

#All Men=1
table3 <- table2[table2[,14]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table3,  Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table3,  Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001889084
anova(model0,model1)[2,7]
#P-value:0.03703822 Significant

coef(model1)[11]
#0.1290108

plot(Effect(focal.predictors = "ScoreInt",model1))


#All Women = 2
table4 <- table2[table2[,14]==2,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table4, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table4, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.00605965
anova(model0,model1)[2,7]
#P-value:0.0001574363 Significant

coef(model1)[11]
#0.2096399

plot(Effect(focal.predictors = "ScoreInt",model1))



#CHILDTRAUMA
plot(table3[,18])
summary(table3[,18])

#All ChildTrauma=1
table5 <- table2[table2[,18]==1,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table5, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table5, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001273307
anova(model0,model1)[2,7]
#P-value:0.09181807 Significant

coef(model1)[12]
#0.09546449

plot(Effect(focal.predictors = "ScoreInt",model1))


#All ChildTrauma = 0
table6 <- table2[table2[,18]==0,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table6, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table6, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.006943723
anova(model0,model1)[2,7]
#P-value:4.125629e-05. Significant

coef(model1)[12]
#0.2487037

plot(Effect(focal.predictors = "ScoreInt",model1))



#CANNABIS
plot(table2[,19])
summary(table2[,19])

#All Cannabis=1
table7 <- table2[table2[,19]==1,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table7, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table7, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.003869763
anova(model0,model1)[2,7]
#P-value:0.01257738 Significant

coef(model1)[12]
#0.1641588

plot(Effect(focal.predictors = "ScoreInt",model1))


#All Cannabis = 0
table8 <- table2[table2[,19]==0,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table8,  Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table8, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004917415
anova(model0,model1)[2,7]
#P-value:0.0001235936 Significant

coef(model1)[12]
#0.2124349

plot(Effect(focal.predictors = "ScoreInt",model1))




############## ANALYSIS OF SIGNIFICANT ENVIRONMENTAL VARIABLES IN RELATION TO GENDER AND COUNTRY IN THE EXPLAINED VARIANCE OF THE SCOREINT -----------------------

########### STRATIFICATION -------------------------------

#SEX

#All Men=1
table3 <- table2[table2[,14]==1,]


#If Cannabis=1
plot(table3[,19])
summary(table3[,19])
table9 <- table3[table3[,19]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table9, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table9, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002125815
anova(model0,model1)[2,7]
#P-value:0.1497436 Not Significant

coef(model1)[11]
#0.1315059

plot(Effect(focal.predictors = "ScoreInt",model1))


#If Cannabis=0
table10 <- table3[table3[,19]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table10,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table10,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002510342
anova(model0,model1)[2,7]
#P-value:0.072816 Significant

coef(model1)[11]
#0.1598209

plot(Effect(focal.predictors = "ScoreInt",model1))



#If ChildTrauma=1
plot(table3[,18])
summary(table3[,18])
table11 <- table3[table3[,18]==1,]



model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table11,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table11,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.003364179
anova(model0,model1)[2,7]
#P-value:0.05675112 Significant

coef(model1)[11]
#0.1646422

plot(Effect(focal.predictors = "ScoreInt",model1))


#If ChildTrauma=0
table12 <- table3[table3[,18]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table12,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table12,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.000921973
anova(model0,model1)[2,7]
#P-value:0.2906303 Not Significant

coef(model1)[11]
#0.09466942

plot(Effect(focal.predictors = "ScoreInt",model1))



#All Women=2
table4 <- table2[table2[,14]==2,]
plot(table4[,19])

#If Cannabis=1
summary(table4[,19])
table13 <- table4[table4[,19]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table13, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table13, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.007076503
anova(model0,model1)[2,7]
#P-value:0.03600792 Significant

coef(model1)[11]
#0.2012478

plot(Effect(focal.predictors = "ScoreInt",model1))


#If Cannabis=0
table14 <- table4[table4[,19]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table14,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table14,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.00669537
anova(model0,model1)[2,7]
#P-value:0.0007390783 Significant

coef(model1)[11]
#0.2408343

plot(Effect(focal.predictors = "ScoreInt",model1))



#If ChildTrauma=1
summary(table4[,18])
table15 <- table4[table4[,18]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table15,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table15,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0004856687
anova(model0,model1)[2,7]
#P-value:0.4572922 Not Significant

coef(model1)[11]
#0.05614165

plot(Effect(focal.predictors = "ScoreInt",model1))


#If ChildTrauma=0
table16 <- table4[table4[,18]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table16,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table16,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01721333
anova(model0,model1)[2,7]
#P-value:6.091931e-06 Significant

coef(model1)[11]
#0.3758889

plot(Effect(focal.predictors = "ScoreInt",model1))



#COUNTRY

#If Men and Cannabis=1
table9 <- table3[table3[,19]==1,]
plot(table9[,16]) #Not enough data
summary(table9[,16])


#If Men and Cannabis=0
table10 <- table3[table3[,19]==0,]
plot(table10[,16]) #Lets analyze only Turkey, Brazil and Spain
summary(table10[,16])

#Turkey
table17 <- table10[table10[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table17,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table17,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.00866487
anova(model0,model1)[2,7]
#P-value:0.0228129 Significant

coef(model1)[11]
#0.8796696

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table18 <- table10[table10[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table18,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table18,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001224935
anova(model0,model1)[2,7]
#P-value:0.6392428 Not Significant

coef(model1)[11]
#0.1166989

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table19 <- table10[table10[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table19,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table19,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.006867442
anova(model0,model1)[2,7]
#P-value:0.377029 Not Significant

coef(model1)[11]
#0.9529289

plot(Effect(focal.predictors = "ScoreInt",model1))



#If Men and ChildTrauma=1
table11 <- table3[table3[,18]==1,]
plot(table11[,16]) #Turkey, Brazil and Spain
summary(table11[,16])

#Turkey
table20 <- table11[table11[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table20,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table20,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01213173
anova(model0,model1)[2,7]
#P-value:0.03949895 Significant

coef(model1)[11]
#1.100256

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table21 <- table11[table11[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table21,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table21,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0009782877
anova(model0,model1)[2,7]
#P-value:0.7893646 Not Significant

coef(model1)[11]
#-0.1066479

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table22 <- table11[table11[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table22,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table22,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002155207
anova(model0,model1)[2,7]
#P-value:0.5620878 Not Significant

coef(model1)[11]
#0.4566211

plot(Effect(focal.predictors = "ScoreInt",model1))


#If Men and ChildTrauma=0
table12 <- table3[table3[,18]==0,]
plot(table12[,16]) #Turkey, Brazil and Spain
summary(table12[,16])

#Turkey
table23 <- table12[table12[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table23,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table23,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002277587
anova(model0,model1)[2,7]
#P-value:0.4273121 Not Significant

coef(model1)[11]
#0.4319485

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table24 <- table12[table12[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table24,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table24,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01573231
anova(model0,model1)[2,7]
#P-value:0.1087101 Significant

coef(model1)[11]
#0.459802

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table25 <- table12[table12[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table25,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table25,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009551057
anova(model0,model1)[2,7]
#P-value:0.06068886 Significant

coef(model1)[11]
#0.8365687

plot(Effect(focal.predictors = "ScoreInt",model1))




#If Women and Cannabis=1
table13 <- table4[table4[,19]==1,]
plot(table13[,16]) #Not enough data
summary(table13[,16])


#If Women and Cannabis=0
table14 <- table4[table4[,19]==0,]
plot(table14[,16])
summary(table14[,16])


#Turkey
table26 <- table14[table14[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table26,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table26,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01089173
anova(model0,model1)[2,7]
#P-value:0.004202748 Not Significant

coef(model1)[11]
#0.9854173

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table27 <- table14[table14[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table27,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table27,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.05566559
anova(model0,model1)[2,7]
#P-value:0.0002558629 Significant

coef(model1)[11]
#0.9461085

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table28 <- table14[table14[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table28,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table28,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.03727711
anova(model0,model1)[2,7]
#P-value:0.01507135 Significant

coef(model1)[11]
#1.826918

plot(Effect(focal.predictors = "ScoreInt",model1))



#If Women and ChildTrauma=1
table15 <- table4[table4[,18]==1,]
plot(table15[,16])
summary(table15[,16])

#Turkey
table29 <- table15[table15[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table29,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table29,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01053855
anova(model0,model1)[2,7]
#P-value:0.03356807 Significant

coef(model1)[11]
#0.9624394

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table30 <- table15[table15[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table30,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table30,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01487376
anova(model0,model1)[2,7]
#P-value:0.2071318 Not Significant

coef(model1)[11]
#0.5818323

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table31 <- table15[table15[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table31,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table31,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.04825082
anova(model0,model1)[2,7]
#P-value:0.02248161 Significant

coef(model1)[11]
#2.066882

plot(Effect(focal.predictors = "ScoreInt",model1))





#If Women and ChildTrauma=0
table16 <- table4[table4[,18]==0,]
plot(table16[,16])
summary(table16[,16])

#Turkey
table32 <- table16[table16[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table32,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table32,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.02307467
anova(model0,model1)[2,7]
#P-value:0.007861864 Significant

coef(model1)[11]
#1.515588

plot(Effect(focal.predictors = "ScoreInt",model1))


#Brazil
table33 <- table16[table16[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table33,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table33,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.09760011
anova(model0,model1)[2,7]
#P-value:5.968723e-05 Significant

coef(model1)[11]
#1.143132

plot(Effect(focal.predictors = "ScoreInt",model1))


#Spain
table34 <- table16[table16[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table34,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreInt, data=table34,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.04569427
anova(model0,model1)[2,7]
#P-value:0.001558229 Not Significant

coef(model1)[11]
#2.066895

plot(Effect(focal.predictors = "ScoreInt",model1))




#############################  STEPWISE BEST MODEL SELECTION -----------------------------------------------------------

#Models with only genetic variables

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(EduLevel ~ Age+Sex+ScoreInt+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data = table2, Hess=T, method='logistic')
model0 = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
model0 <- polr(EduLevel ~ Age + Sex + ScoreInt + PCA2 + PCA4 + 
                 PCA5 + PCA9 + Age:PCA2 + ScoreInt:PCA5 + ScoreInt:PCA4 + 
                 Sex:PCA2 + Age:Sex + Age:ScoreInt + PCA2:PCA4 + PCA2:PCA9, data = table2, Hess=T, method='logistic')
extractAIC(model0)[2] 
#4573.865
var1 <- 1-(model0$deviance/model_null$deviance)
var1
#R2:0.06011746

coef(model0)[3]
#0.5844347


#Models with environmental factors 

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(EduLevel ~ Age+Sex+ScoreInt+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+Country+ChildTrauma+Cannabis, data = table2, Hess=T, method='logistic')
best_model = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
best_model <- polr(EduLevel ~ Age + ScoreInt + PCA5 + PCA7 + PCA8 + 
                     PCA10 + Country + ChildTrauma + Cannabis + Age:Country + 
                     ScoreInt:Country + Age:Cannabis + PCA7:PCA8 + PCA5:Cannabis + 
                     Age:PCA7 + Age:PCA10, data = table2, Hess=T, method='logistic')
extractAIC(best_model)[2] 
#4475.786
var1 <- 1-(best_model$deviance/model_null$deviance)
var1
#R2:0.08704763

coef(best_model)[2]
#1.108015



####################################### COGNITION SCORES ---------------------------------------------------------------------


########### ANALISIS OF EACH ENVIRONMENTAL VARIABLE IN THE EFFECT OF SCOREINT ------------------------------


########### STRATIFICATION -------------------------------

#SEX
plot(table2[,14])

#All Men=1
table3 <- table2[table2[,14]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table3,  Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table3,  Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001941437
anova(model0,model1)[2,7]
#P-value:0.03451178 Significant

coef(model1)[11]
#0.1227792

plot(Effect(focal.predictors = "ScoreCog",model1))


#All Women = 2
table4 <- table2[table2[,14]==2,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table4, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table4, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001313686
anova(model0,model1)[2,7]
#P-value:0.07848315 Significant

coef(model1)[11]
#0.1037948

plot(Effect(focal.predictors = "ScoreCog",model1))



#CHILDTRAUMA
plot(table3[,18])

#All ChildTrauma=1
table5 <- table2[table2[,18]==1,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table5, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table5, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001546111
anova(model0,model1)[2,7]
#P-value:0.06320806 Significant

coef(model1)[12]
#0.1061947

plot(Effect(focal.predictors = "ScoreCog",model1))


#All ChildTrauma = 0
table6 <- table2[table2[,18]==0,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table6, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table6, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002422128
anova(model0,model1)[2,7]
#P-value:0.01544793 Significant

coef(model1)[12]
#0.1462727

plot(Effect(focal.predictors = "ScoreCog",model1))



#CANNABIS
plot(table2[,19])

#All Cannabis=1
table7 <- table2[table2[,19]==1,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table7, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table7, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004101391
anova(model0,model1)[2,7]
#P-value:0.0101958 Not significant

coef(model1)[12]
#0.1770802

plot(Effect(focal.predictors = "ScoreCog",model1))


#All Cannabis = 0
table8 <- table2[table2[,19]==0,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table8,  Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table8, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0004708191
anova(model0,model1)[2,7]
#P-value:0.2348901 Not Significant

coef(model1)[12]
#0.06333791

plot(Effect(focal.predictors = "ScoreCog",model1))




############## ANALYSIS OF SIGNIFICANT ENVIRONMENTAL VARIABLES IN RELATION TO GENDER AND COUNTRY IN THE EXPLAINED VARIANCE OF THE ScoreCog -----------------------

########### STRATIFICATION -------------------------------

#SEX

#All Men=1
table3 <- table2[table2[,14]==1,]
plot(table3[,19])

#If Cannabis=1
table9 <- table3[table3[,19]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table9, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table9, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.007177526
anova(model0,model1)[2,7]
#P-value:0.008126029 Significant

coef(model1)[11]
#0.230781

plot(Effect(focal.predictors = "ScoreCog",model1))


#If Cannabis=0
table10 <- table3[table3[,19]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table10,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table10,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0001871868
anova(model0,model1)[2,7]
#P-value:0.6242192 Not Significant

coef(model1)[11]
#0.03950657

plot(Effect(focal.predictors = "ScoreCog",model1))



#If ChildTrauma=1
table11 <- table3[table3[,18]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table11,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table11,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001213728
anova(model0,model1)[2,7]
#P-value:0.2524696 Not Significant

coef(model1)[11]
#0.09119285

plot(Effect(focal.predictors = "ScoreCog",model1))


#If ChildTrauma=0
table12 <- table3[table3[,18]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table12,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table12,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.003511449
anova(model0,model1)[2,7]
#P-value:0.03917903 Not Significant

coef(model1)[11]
#0.1757464

plot(Effect(focal.predictors = "ScoreCog",model1))



#All Women=2
table4 <- table2[table2[,14]==2,]
plot(table4[,19])

#If Cannabis=1
table13 <- table4[table4[,19]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table13, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table13, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0008914769
anova(model0,model1)[2,7]
#P-value:0.4567338 Not Significant

coef(model1)[11]
#0.08610457

plot(Effect(focal.predictors = "ScoreCog",model1))


#If Cannabis=0
table14 <- table4[table4[,19]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table14,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table14,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0007729787
anova(model0,model1)[2,7]
#P-value:0.2515316 Not Significant

coef(model1)[11]
#0.08178457

plot(Effect(focal.predictors = "ScoreCog",model1))



#If ChildTrauma=1
table15 <- table4[table4[,18]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table15,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table15,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001616809
anova(model0,model1)[2,7]
#P-value:0.1750284 Not Significant

coef(model1)[11]
#0.1116297

plot(Effect(focal.predictors = "ScoreCog",model1))


#If ChildTrauma=0
table16 <- table4[table4[,18]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table16,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table16,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001684129
anova(model0,model1)[2,7]
#P-value:0.1571238 Not Significant

coef(model1)[11]
#0.1222054

plot(Effect(focal.predictors = "ScoreCog",model1))



#COUNTRY

#If Men and Cannabis=1
table9 <- table3[table3[,19]==1,]
plot(table9[,16]) #Not enough data



#If Men and Cannabis=0
table10 <- table3[table3[,19]==0,]
plot(table10[,16]) #Lets analyze only Turkey, Brazil and Spain

#Turkey
table17 <- table10[table10[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table17,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table17,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002868364
anova(model0,model1)[2,7]
#P-value:0.1902543 Not Significant

coef(model1)[11]
#0.2179558

plot(Effect(focal.predictors = "ScoreCog",model1))


#Brazil
table18 <- table10[table10[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table18,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table18,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0009850493
anova(model0,model1)[2,7]
#P-value:0.6742227 Not Significant

coef(model1)[11]
#-0.07049409 

plot(Effect(focal.predictors = "ScoreCog",model1))


#Spain
table19 <- table10[table10[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table19,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table19,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.02544945
anova(model0,model1)[2,7]
#P-value:0.08902623 Not Significant

coef(model1)[11]
#0.6715829

plot(Effect(focal.predictors = "ScoreCog",model1))



#If Men and ChildTrauma=1
table11 <- table3[table3[,18]==1,]
plot(table11[,16]) #Turkey, Brazil and Spain


#Turkey
table20 <- table11[table11[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table20,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table20,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009963434
anova(model0,model1)[2,7]
#P-value:0.06205485 Significant

coef(model1)[11]
#0.4075072

plot(Effect(focal.predictors = "ScoreCog",model1))


#Brazil
table21 <- table11[table11[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table21,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table21,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001037574
anova(model0,model1)[2,7]
#P-value:0.7832308 Not Significant

coef(model1)[11]
#-0.06989028

plot(Effect(focal.predictors = "ScoreCog",model1))


#Spain
table22 <- table11[table11[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table22,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table22,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004317577
anova(model0,model1)[2,7]
#P-value:0.411896 Not Significant

coef(model1)[11]
#0.2494714

plot(Effect(focal.predictors = "ScoreCog",model1))


#If Men and ChildTrauma=0
table12 <- table3[table3[,18]==0,]
plot(table12[,16]) #Turkey, Brazil and Spain


#Turkey
table23 <- table12[table12[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table23,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table23,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001176353
anova(model0,model1)[2,7]
#P-value:0.5683506 Not Significant

coef(model1)[11]
#-0.1444779

plot(Effect(focal.predictors = "ScoreCog",model1))


#Brazil
table24 <- table12[table12[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table24,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table24,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0002190568
anova(model0,model1)[2,7]
#P-value:0.849878 Not Significant

coef(model1)[11]
#-0.04109466

plot(Effect(focal.predictors = "ScoreCog",model1))


#Spain
table25 <- table12[table12[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table25,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table25,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009473991
anova(model0,model1)[2,7]
#P-value:0.06173802 Significant

coef(model1)[11]
#0.3732841

plot(Effect(focal.predictors = "ScoreCog",model1))




#If Women and Cannabis=1
table13 <- table4[table4[,19]==1,]
plot(table13[,16]) #Not enough data



#If Women and Cannabis=0
table14 <- table4[table4[,19]==0,]
plot(table14[,16])


#Turkey
table26 <- table14[table14[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table26,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table26,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002894468
anova(model0,model1)[2,7]
#P-value:0.1400357 Not Significant

coef(model1)[11]
#0.243506

plot(Effect(focal.predictors = "ScoreCog",model1))


#Brazil
table27 <- table14[table14[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table27,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table27,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0001010916
anova(model0,model1)[2,7]
#P-value:0.8761791 Significant

coef(model1)[11]
#-0.0270017

plot(Effect(focal.predictors = "ScoreCog",model1))


#Spain
table28 <- table14[table14[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table28,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table28,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0192708
anova(model0,model1)[2,7]
#P-value:0.08052574 Significant

coef(model1)[11]
#0.5196619

plot(Effect(focal.predictors = "ScoreCog",model1))



#If Women and ChildTrauma=1
table15 <- table4[table4[,18]==1,]
plot(table16[,16])


#Turkey
table29 <- table15[table15[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table29,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table29,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001868215
anova(model0,model1)[2,7]
#P-value:0.3708924 Not Significant

coef(model1)[11]
#0.1985456

plot(Effect(focal.predictors = "ScoreCog",model1))


#Brazil
table30 <- table15[table15[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table30,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table30,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001018122
anova(model0,model1)[2,7]
#P-value:0.7413655 Not Significant

coef(model1)[11]
#0.08396447

plot(Effect(focal.predictors = "ScoreCog",model1))


#Spain
table31 <- table15[table15[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table31,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table31,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0154765
anova(model0,model1)[2,7]
#P-value:0.1961902 Not Significant

coef(model1)[11]
#0.5681977

plot(Effect(focal.predictors = "ScoreCog",model1))





#If Women and ChildTrauma=0
table16 <- table4[table4[,18]==0,]
plot(table16[,16])


#Turkey
table32 <- table16[table16[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table32,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table32,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.006409069
anova(model0,model1)[2,7]
#P-value:0.1612742 Significant

coef(model1)[11]
#0.359629

plot(Effect(focal.predictors = "ScoreCog",model1))


#Brazil
table33 <- table16[table16[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table33,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table33,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.02606932
anova(model0,model1)[2,7]
#P-value:0.03802917 Significant

coef(model1)[11]
#-0.4601937

plot(Effect(focal.predictors = "ScoreCog",model1))


#Spain
table34 <- table16[table16[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table34,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreCog, data=table34,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.007265829
anova(model0,model1)[2,7]
#P-value:0.2071204 Not Significant

coef(model1)[11]
#0.3129403

plot(Effect(focal.predictors = "ScoreCog",model1))




#############################  STEPWISE BEST MODEL SELECTION -----------------------------------------------------------

#Models with only genetic variables

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(EduLevel ~ Age+Sex+ScoreCog+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data = table2, Hess=T, method='logistic')
model0 = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
model0 <- polr(EduLevel ~ Age + ScoreCog + PCA2 + PCA4 + PCA5 + 
                 PCA9 + Age:PCA2 + PCA2:PCA4 + PCA2:PCA9 + Age:ScoreCog, data = table2, Hess=T, method='logistic')
extractAIC(model0)[2] 
#4623.466
var0 <- 1-(model0$deviance/model_null$deviance)
var0
#0.04777833

coef(model0)[2]
#-0.2181475

#Models with environmental factors 

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(EduLevel ~ Age+Sex+ScoreCog+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+Country+ChildTrauma+Cannabis, data = table2, Hess=T, method='logistic')
best_model = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
best_model <- polr(EduLevel ~ Age + Sex + ScoreCog + PCA2 + PCA3 + 
                     PCA5 + PCA7 + PCA8 + PCA10 + Country + ChildTrauma + Cannabis + 
                     Age:Country + PCA5:Country + Age:Cannabis + Sex:PCA2 + Age:Sex + 
                     Age:PCA7 + Age:PCA3 + ScoreCog:PCA8 + PCA2:PCA10 + PCA5:Cannabis, data = table2, Hess=T, method='logistic')
extractAIC(best_model)[2] 
#4486.534
var0 <- 1-(best_model$deviance/model_null$deviance)
var0
#0.08730674

coef(best_model)[3]
#0.1083015




####################################### EDUCATIONAL ATTAINMENT SCORES ---------------------------------------------------------------------


########### ANALISIS OF EACH ENVIRONMENTAL VARIABLE IN THE EFFECT OF SCOREEA ------------------------------


########### STRATIFICATION -------------------------------

#SEX
plot(table2[,14])

#All Men=1
table3 <- table2[table2[,14]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table3,  Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table3,  Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004740347
anova(model0,model1)[2,7]
#P-value:0.0009553129 Not significant

coef(model1)[11]
#0.1977722

plot(Effect(focal.predictors = "ScoreEA",model1))


#All Women = 2
table4 <- table2[table2[,14]==2,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table4, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table4, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009441788
anova(model0,model1)[2,7]
#P-value:2.391106e-06 Significant

coef(model1)[11]
#0.2704709

plot(Effect(focal.predictors = "ScoreEA",model1))



#CHILDTRAUMA
plot(table3[,18])

#All ChildTrauma=1
table5 <- table2[table2[,18]==1,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table5, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table5, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.004014686
anova(model0,model1)[2,7]
#P-value:0.002757496 Significant

coef(model1)[12]
#0.1707827

plot(Effect(focal.predictors = "ScoreEA",model1))


#All ChildTrauma = 0
table6 <- table2[table2[,18]==0,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table6, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table6, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.009317442
anova(model0,model1)[2,7]
#P-value:2.036669e-06 Significant

coef(model1)[12]
#0.2863316

plot(Effect(focal.predictors = "ScoreEA",model1))



#CANNABIS
plot(table2[,19])

#All Cannabis=1
table7 <- table2[table2[,19]==1,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table7, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table7, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.008762669
anova(model0,model1)[2,7]
#P-value:0.0001731825 Significant

coef(model1)[12]
#0.2486167

plot(Effect(focal.predictors = "ScoreEA",model1))


#All Cannabis = 0
table8 <- table2[table2[,19]==0,]

model0 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table8,  Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+Sex+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table8, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.007169615
anova(model0,model1)[2,7]
#P-value:3.562933e-06 Significant

coef(model1)[12]
#0.2547889

plot(Effect(focal.predictors = "ScoreEA",model1))




############## ANALYSIS OF SIGNIFICANT ENVIRONMENTAL VARIABLES IN RELATION TO GENDER AND COUNTRY IN THE EXPLAINED VARIANCE OF THE ScoreEA -----------------------

########### STRATIFICATION -------------------------------


#SEX

#All Men=1
table3 <- table2[table2[,14]==1,]
plot(table3[,19])

#If Cannabis=1
table9 <- table3[table3[,19]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table9, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table9, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.00631554
anova(model0,model1)[2,7]
#P-value:0.01303641 Significant

coef(model1)[11]
#0.2250243

plot(Effect(focal.predictors = "ScoreEA",model1))


#If Cannabis=0
table10 <- table3[table3[,19]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table10,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table10,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.005290033
anova(model0,model1)[2,7]
#P-value:0.009207842 Significant

coef(model1)[11]
#0.2187674

plot(Effect(focal.predictors = "ScoreEA",model1))



#If ChildTrauma=1
table11 <- table3[table3[,18]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table11,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table11,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.006699677
anova(model0,model1)[2,7]
#P-value:0.007174309 Significant

coef(model1)[11]
#0.2251275

plot(Effect(focal.predictors = "ScoreEA",model1))


#If ChildTrauma=0
table12 <- table3[table3[,18]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table12,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table12,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.003033993
anova(model0,model1)[2,7]
#P-value:0.05524118 Not Significant

coef(model1)[11]
#0.1656231

plot(Effect(focal.predictors = "ScoreEA",model1))



#All Women=2
table4 <- table2[table2[,14]==2,]
plot(table4[,19])

#If Cannabis=1
table13 <- table4[table4[,19]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table13, Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table13, Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0136984
anova(model0,model1)[2,7]
#P-value:0.003530045 Significant

coef(model1)[11]
#0.2866469

plot(Effect(focal.predictors = "ScoreEA",model1))


#If Cannabis=0
table14 <- table4[table4[,19]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table14,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table14,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.008661559
anova(model0,model1)[2,7]
#P-value:0.0001238827 Significant

coef(model1)[11]
#0.2831879

plot(Effect(focal.predictors = "ScoreEA",model1))



#If ChildTrauma=1
table15 <- table4[table4[,18]==1,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table15,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table15,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002881683
anova(model0,model1)[2,7]
#P-value:0.07020135 Significant

coef(model1)[11]
#0.1421121

plot(Effect(focal.predictors = "ScoreEA",model1))


#If ChildTrauma=0
table16 <- table4[table4[,18]==0,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table16,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table16,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01955559
anova(model0,model1)[2,7]
#P-value:1.427678e-06 Significant

coef(model1)[11]
#0.4130979

plot(Effect(focal.predictors = "ScoreEA",model1))



#COUNTRY

#If Men and Cannabis=1
table9 <- table3[table3[,19]==1,]
plot(table9[,16]) #Not enough data



#If Men and Cannabis=0
table10 <- table3[table3[,19]==0,]
plot(table10[,16]) #Lets analyze only Turkey, Brazil and Spain

#Turkey
table17 <- table10[table10[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table17,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table17,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.00492078
anova(model0,model1)[2,7]
#P-value:0.08623689 Significant

coef(model1)[11]
#0.4936515

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table18 <- table10[table10[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table18,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table18,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001002923
anova(model0,model1)[2,7]
#P-value:0.6714519 Not Significant

coef(model1)[11]
#0.09455561

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table19 <- table10[table10[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table19,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table19,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.03482863
anova(model0,model1)[2,7]
#P-value:0.04665808 Significant

coef(model1)[11]
#1.350775

plot(Effect(focal.predictors = "ScoreEA",model1))



#If Men and ChildTrauma=1
table11 <- table3[table3[,18]==1,]
plot(table11[,16]) #Turkey, Brazil and Spain


#Turkey
table20 <- table11[table11[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table20,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table20,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0004412629
anova(model0,model1)[2,7]
#P-value:0.6945594 Not Significant

coef(model1)[11]
#0.1436418

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table21 <- table11[table11[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table21,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table21,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0117871
anova(model0,model1)[2,7]
#P-value:0.3537907 Not Significant

coef(model1)[11]
#0.3751365

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table22 <- table11[table11[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table22,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table22,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01313695
anova(model0,model1)[2,7]
#P-value:0.1523369 Not Significant

coef(model1)[11]
#0.8027122

plot(Effect(focal.predictors = "ScoreEA",model1))


#If Men and ChildTrauma=0
table12 <- table3[table3[,18]==0,]
plot(table12[,16]) #Turkey, Brazil and Spain


#Turkey
table23 <- table12[table12[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table23,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table23,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.008623562
anova(model0,model1)[2,7]
#P-value:0.1224424 Not Significant

coef(model1)[11]
#0.7412459

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table24 <- table12[table12[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table24,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table24,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002425155
anova(model0,model1)[2,7]
#P-value:0.5288444 Not Significant

coef(model1)[11]
#0.1603829

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table25 <- table12[table12[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table25,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table25,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.01617141
anova(model0,model1)[2,7]
#P-value:0.01465655 Significant

coef(model1)[11]
#0.6802892

plot(Effect(focal.predictors = "ScoreEA",model1))




#If Women and Cannabis=1
table13 <- table4[table4[,19]==1,]
plot(table13[,16]) #Not enough data



#If Women and Cannabis=0
table14 <- table4[table4[,19]==0,]
plot(table14[,16])


#Turkey
table26 <- table14[table14[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table26,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table26,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.001092978
anova(model0,model1)[2,7]
#P-value:0.3645174 Not Significant

coef(model1)[11]
#0.2507565

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table27 <- table14[table14[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table27,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table27,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0514151
anova(model0,model1)[2,7]
#P-value:0.0004414897 Significant

coef(model1)[11]
#0.8459277

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table28 <- table14[table14[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table28,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table28,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.08877218
anova(model0,model1)[2,7]
#P-value:0.0001761644 Significant

coef(model1)[11]
#2.005725

plot(Effect(focal.predictors = "ScoreEA",model1))



#If Women and ChildTrauma=1
table15 <- table4[table4[,18]==1,]
plot(table16[,16])


#Turkey
table29 <- table15[table15[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table29,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table29,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.002882462
anova(model0,model1)[2,7]
#P-value:0.2663683 Not Significant

coef(model1)[11]
#0.4270975

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table30 <- table15[table15[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table30,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table30,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.04147696
anova(model0,model1)[2,7]
#P-value:0.03515422 Significant

coef(model1)[11]
#0.9277485

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table31 <- table15[table15[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table31,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table31,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.08365195
anova(model0,model1)[2,7]
#P-value:0.00265691 Significant

coef(model1)[11]
#2.358086

plot(Effect(focal.predictors = "ScoreEA",model1))





#If Women and ChildTrauma=0
table16 <- table4[table4[,18]==0,]
plot(table16[,16])


#Turkey
table32 <- table16[table16[,16]==10,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table32,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table32,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.0001897492
anova(model0,model1)[2,7]
#P-value:0.8095331 Not Significant

coef(model1)[11]
#0.1099066

plot(Effect(focal.predictors = "ScoreEA",model1))


#Brazil
table33 <- table16[table16[,16]==12,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table33,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table33,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.07483244
anova(model0,model1)[2,7]
#P-value:0.0004400678 Significant

coef(model1)[11]
#0.9764185

plot(Effect(focal.predictors = "ScoreEA",model1))


#Spain
table34 <- table16[table16[,16]==7,]

model0 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data=table34,Hess=T, method='logistic')
model1 <- polr(EduLevel~Age+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+ScoreEA, data=table34,Hess=T, method='logistic')
summary(model0)
summary(model1)

var1 <- 1-(model1$deviance/model0$deviance)
var1
#R2:0.06330188
anova(model0,model1)[2,7]
#P-value:0.0001964155 Significant

coef(model1)[11]
#1.818471

plot(Effect(focal.predictors = "ScoreEA",model1))




#############################  STEPWISE BEST MODEL SELECTION -----------------------------------------------------------

#Models with only genetic variables

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(EduLevel ~ Age+Sex+ScoreEA+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10, data = table2, Hess=T, method='logistic')
model0 = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
model0 <- polr(EduLevel ~ Age + Sex + ScoreEA + PCA2 + PCA3 + 
                 PCA4 + PCA5 + PCA8 + Age:PCA2 + ScoreEA:PCA3 + Sex:PCA2 + 
                 Age:Sex + PCA2:PCA3 + Sex:PCA5, data = table2)
extractAIC(model0)[2] 
#4586.261
var0 <- 1-(model0$deviance/model_null$deviance)
var0
#0.05713713

coef(model0)[3]
#0.3699367


#Models with environmental factors 

# selection with AIC criterion:
set.seed(3890)
fit.all <- polr(EduLevel ~ Age+Sex+ScoreEA+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+Country+ChildTrauma+Cannabis, data = table2, Hess=T, method='logistic')
best_model = stepAIC(fit.all, scope =  . ~ .^2, trace= F, steps = 1000, k=qchisq(0.05, 1, lower.tail = F)) # is the fit model for this data considering interactions
best_model <- polr(EduLevel ~ Age + Sex + ScoreEA + PCA2 + PCA3 + 
                     PCA7 + PCA8 + PCA10 + Country + ChildTrauma + Cannabis + 
                     Age:Country + ScoreEA:Country + Age:Cannabis + PCA2:PCA3 + 
                     Sex:PCA2 + Age:Sex + Age:PCA7 + PCA7:PCA8 + Age:PCA10, data = table2, Hess=T, method='logistic')
extractAIC(best_model)[2] 
#4449.202
var0 <- 1-(best_model$deviance/model_null$deviance)
var0
#0.09420745

coef(best_model)[3]
#0.4159375

#Best model is with ScoreEA



#########################  PLOTING THE EFFECTS OF THE VARIABLES IN THE BEST MODEL ------------------------------------------------------

plot(Effect(focal.predictors = "Age",best_model))
#Older people have more probability to have higher EduLevels

plot(Effect(focal.predictors = "Sex",best_model))
#The women have much lower probability to have higher EduLevels

plot(Effect(focal.predictors = "Country",best_model))

plot(Effect(focal.predictors = "ScoreEA",best_model))
#The higher the Score, the higher the EduLevels

plot(Effect(focal.predictors = "ChildTrauma",best_model))
#The presence of ChildTrauma involves to have less probability to gain higher EduLevels

plot(Effect(focal.predictors = "Cannabis",best_model))
#The consumption of Cannabis involves to have a little bit more probability to gain higher EduLevels



################# PREDICTION POWER OF THE BEST MODEL -------------------------------------------

#SENSIBILITY TEST
ROC(form=EduLevel ~ Age + Sex + ScoreEA + PCA2 + PCA3 + 
      PCA7 + PCA8 + PCA10 + Country + ChildTrauma + Cannabis + 
      Age:Country + ScoreEA:Country + Age:Cannabis + PCA2:PCA3 + 
      Sex:PCA2 + Age:Sex + Age:PCA7 + PCA7:PCA8 + Age:PCA10, data=table2, plot="ROC", lwd=3, cex=1.5, AUC=TRUE)
#AUC=0.75
#Acceptable AUC for classification with 3 levels



#CONFUSION MATRIX

#Prepare Training and Test Data
set.seed(3890)
trainingRows <- sample(1:nrow(table2), 0.7 * nrow(table2))
trainingData <- table2[trainingRows, ]
testData <- table2[-trainingRows, ]

finalmod2 <- polr(EduLevel ~ Age + Sex + ScoreEA + PCA2 + PCA3 + 
                    PCA7 + PCA8 + PCA10 + Country + ChildTrauma + Cannabis + 
                    Age:Country + ScoreEA:Country + Age:Cannabis + PCA2:PCA3 + 
                    Sex:PCA2 + Age:Sex + Age:PCA7 + PCA7:PCA8 + Age:PCA10, data=table2, Hess=TRUE, method="logistic")
summary(finalmod2)

predictedPhenotype <- predict(finalmod2, testData)  # predict the classes directly
head(predictedPhenotype)

predictedScores <- predict(finalmod2, testData, type="p")  # predict the probabilites
head(predictedScores)

ConfMat <- table(testData$EduLevel, predictedPhenotype)  # confusion matrix
ConfMat

mean(as.character(testData$EduLevel) != as.character(predictedPhenotype))  # missclassification error
#A 47% of error.


n = length(testData$EduLevel)
accuracy <- sum(diag(ConfMat)) / n
accuracy
#A 53% of accuracy.



