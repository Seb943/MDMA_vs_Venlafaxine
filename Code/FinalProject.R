###################Data Analysis and Exploration Project###################

#Sébastien CARARO 2019 - University of Trento

#The dataset, which provides results from testing both MDMA and VLX on Dark Agouti rats, 
#we could like to question whether or not certain zone of the body are more affected by such
#treatments (3 zones were tested : Frontal Cortex, Hypocampus, Dorsal Raphe). We could also 
#reflect on whether the effects of drugs are cumulative. 

cat("\014") 

list.files()

rm(list=ls())

ls()

getwd()

setwd(dir = "C:/Users/Sébastien CARARO/Documents")

#??? grid.arrange(p1,p2,p3) to display multiple plots at once

# Ctrl + Shift + C met la selection en commentaires

#putting the labels on a plot 
#text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.5, labels = colnames())

#Save a data.frame file 
#write.table(data, "data.csv", row.names=FALSE, sep="t",dec=",", na=" ")


#Installing packages
# install.packages("BiocManager")
# BiocManager::install("GEOquery")
# BiocManager::install("useful")
# BiocManager::install("ALL"); data(ALL)
# 
# install.packages("dplyr")


#Include libraries
library("dplyr")

library("GEOquery")
library("useful")
library("ALL")
library("randomForest")
library("caret")
library("genefilter")

#  Loading the file, clear the bad cells, set the data as expressionset and displays the columns names, 
#  then rename then in order to make them explicitly interpretable
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE, parseCharacteristics = FALSE)
length(gse)
show(gse)
na.omit(gse)
#gse <- gse[,-c(1,19,23,24,41)]
gse<-new("ExpressionSet", exprs=as.matrix(gse)) #replaces the exprs() function
head(exprs(gse))
length(exprs(gse)) 
ex <- exprs(gse)
dim(ex)
colnames(ex)
#ex[1:5,]

#Loading the description file to get the Illumina titles
descrip <- read.table(file = "GEO Accession viewer_GSE47541_FULL_modified2.txt", header = TRUE, fill = TRUE)
colnames(descrip)
head(descrip)
dim(descrip)
SYMBOL <- as.factor(descrip[,3])

#Get the signatures
signa <- read.table(file = "SCUDO_testing_signatures(overalldata-ColumnsNotWellAssigned).csv", sep=";")
dim(signa)
SIGNA <- as.factor(signa[c(-1,-2),c(17,18)])
head(SIGNA)
dim(SIGNA)

#Setting the columns names considering the treatment and the zone tested 
colnames(ex) <- c("Fc_MDMA/VLX_1","Fc_MDMA/VLX_2","Fc_MDMA/VLX_3","Fc_MDMA/SAL_1","Fc_MDMA/SAL_2",
                  "Fc_MDMA/SAL_3","Fc_MDMA/SAL_4","Fc_SAL/VLX_1","Fc_SAL/VLX_2","Fc_SAL/VLX_3",
                  "Fc_SAL/VLX_4","Fc_SAL/SAL_1","Fc_SAL/SAL_2","Fc_SAL/SAL_3","Fc_SAL/SAL_4",
                  "Hp_MDMA/VLX_1","Hp_MDMA/VLX_2","Hp_MDMA/VLX_3","Hp_MDMA/VLX_4","Hp_MDMA/SAL_1",
                  "Hp_MDMA/SAL_2","Hp_MDMA/SAL_3","Hp_MDMA/SAL_4","Hp_SAL/VLX_1","Hp_SAL/VLX_2",
                  "Hp_SAL/VLX_3","Hp_SAL/VLX_4","Hp_SAL/SAL_1","Hp_SAL/SAL_2","Hp_SAL/SAL_3",
                  "Hp_SAL/SAL_4","Dr_MDMA/VLX_1","Dr_MDMA/VLX_2","Dr_MDMA/VLX_3","Dr_MDMA/VLX_4",
                  "Dr_MDMA/SAL_1","Dr_MDMA/SAL_3","Dr_MDMA/SAL_4","Dr_SAL/VLX_1","Dr_SAL/VLX_2",
                  "Dr_SAL/VLX_3","Dr_SAL/VLX_4","Dr_SAL/SAL_1","Dr_SAL/SAL_2","Dr_SAL/SAL_3","Dr_SAL/SAL_4")
colnames(ex)
str(ex)

colnames(gse) <- c("Fc_MDMA/VLX_1","Fc_MDMA/VLX_2","Fc_MDMA/VLX_3","Fc_MDMA/SAL_1","Fc_MDMA/SAL_2",
                  "Fc_MDMA/SAL_3","Fc_MDMA/SAL_4","Fc_SAL/VLX_1","Fc_SAL/VLX_2","Fc_SAL/VLX_3",
                  "Fc_SAL/VLX_4","Fc_SAL/SAL_1","Fc_SAL/SAL_2","Fc_SAL/SAL_3","Fc_SAL/SAL_4",
                  "Hp_MDMA/VLX_1","Hp_MDMA/VLX_2","Hp_MDMA/VLX_3","Hp_MDMA/VLX_4","Hp_MDMA/SAL_1",
                  "Hp_MDMA/SAL_2","Hp_MDMA/SAL_3","Hp_MDMA/SAL_4","Hp_SAL/VLX_1","Hp_SAL/VLX_2",
                  "Hp_SAL/VLX_3","Hp_SAL/VLX_4","Hp_SAL/SAL_1","Hp_SAL/SAL_2","Hp_SAL/SAL_3",
                  "Hp_SAL/SAL_4","Dr_MDMA/VLX_1","Dr_MDMA/VLX_2","Dr_MDMA/VLX_3","Dr_MDMA/VLX_4",
                  "Dr_MDMA/SAL_1","Dr_MDMA/SAL_3","Dr_MDMA/SAL_4","Dr_SAL/VLX_1","Dr_SAL/VLX_2",
                  "Dr_SAL/VLX_3","Dr_SAL/VLX_4","Dr_SAL/SAL_1","Dr_SAL/SAL_2","Dr_SAL/SAL_3","Dr_SAL/SAL_4")
colnames(gse)
str(gse)

#Splitting into each zone : Frontal Cortex, Hypocampus, Dorsal Raphe
ex_Fc <- ex[,1:15]
gse_Fc <- gse[,1:15]
colnames(ex_Fc)
dim(ex_Fc)
ex_Hp <- ex[,16:31]
gse_Hp <- gse[,16:31]
colnames(ex_Hp)
dim(ex_Hp)
ex_Dr <- ex[,32:46]
gse_Dr <- gse[,32:46]
colnames(ex_Dr)
dim(ex_Dr)

#Splitting into treatments : MDMA or SALINE, SALINE or VLX
ex_MDMAVLX <- ex[,c(1,2,3,16,17,18,19,32,33,34,35)]
colnames(ex_MDMAVLX)
dim(ex_MDMAVLX)
ex_MDMASAL <- ex[,c(4,5,6,7,20,21,22,23,36,37,38)]
colnames(ex_MDMASAL)
dim(ex_MDMASAL)
ex_SALVLX <- ex[,c(8,9,10,11,24,25,26,27,39,40,41,42)]
colnames(ex_SALVLX)
dim(ex_SALVLX)
ex_SALSAL <- ex[,c(12,13,14,15,28,29,30,31,43,44,45,46)]
colnames(ex_SALSAL)
dim(ex_SALSAL)

#Creating a dataset containing MDMAVLX and SALSAL : 23 samples 
ex_TEST <- ex[,c(1,2,3,16,17,18,19,32,33,34,35,12,13,14,15,28,29,30,31,43,44,45,46)]
gse_TEST <- gse[,c(1,2,3,16,17,18,19,32,33,34,35,12,13,14,15,28,29,30,31,43,44,45,46)]
#we could also remove the outliers : 1,19,23,24,41 : 
#ex_TEST <- ex[,c(2,3,16,17,18,32,33,34,35,12,13,14,15,28,29,30,31,43,44,45,46)]
colnames(ex_TEST)

#Creating a dataset containing SALVLX and SALSAL : 24 samples (12/12)
ex_VLX <- ex[,c(8,9,10,11,24,25,26,27,39,40,41,42,12,13,14,15,28,29,30,31,43,44,45,46)]
gse_VLX <- gse[,c(8,9,10,11,24,25,26,27,39,40,41,42,12,13,14,15,28,29,30,31,43,44,45,46)]
colnames(ex_VLX)
dim(ex_VLX)

#Creating a dataset containing MDMASAL and SALSAL : 23 samples(11/12)
ex_MDMA <- ex[,c(4,5,6,7,20,21,22,23,36,37,38,12,13,14,15,28,29,30,31,43,44,45,46)]
gse_MDMA <- gse[,c(4,5,6,7,20,21,22,23,36,37,38,12,13,14,15,28,29,30,31,43,44,45,46)]
colnames(ex_MDMA)


#removing the outliers on gse and ex
#gse <- gse[,-c(1,19,23,24,41)]
#ex <- ex[,-c(1,19,23,24,41)]


# ex <- ex_TEST
# 
# ex <- ex_Hp
# ex <- ex_Fc
# ex <- ex_Dr
# 
# ex <- ex_MDMAVLX
# ex <- ex_MDMASAL
# ex <- ex_SALVLX
# ex <- ex_SALSAL














############################## WEEK 2 : Data Visualization #################################
#install.packages("BiocManager")
#BiocManager::install("GEOquery")
#library("GEOquery")
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
# getGEO returns a list of expression objects, but ...
length(gse)
# ... shows us there is only one object in it. We assign it to
# the same variable.
#gse <- gse[[1]]
# show what we have:
show(gse)
na.omit(gse)
# The actual expression data are accessible in the "exprs" ... an
# Expression Set, the generic data class that BioConductor uses
# for expression data
gse<-new("ExpressionSet", exprs=as.matrix(gse)) #replaces the exprs() function
head(exprs(gse))
length(exprs(gse)) # 32321 rows times 18 columns
# exprs(gse) is a matrix that we can assign to its own variable, to
# conveniently access the data rows and columns
ex <- exprs(gse)
dim(ex)
colnames(ex)
ex[1:5,]



# #Scale normalization
# # compute MAD, then divide by column
# medians=apply(X,2,median)
# Y=sweep(X,2,medians,"-")
# mad=apply(abs(Y),2,median)
# const=prod(mad)^(1/length(mad))
# scale.normalized.X=sweep(X,2,const/mad,"*")

# simpler approach
scale.normalized.X <- scale(X)

# Analyze value distributions
boxplot(ex)

#Processing raw data : see on slides week 2 for datatypes Affimetrix, Illumina, Agilent
















############################## WEEK 3 : PCA  #################################
#install.packages("BiocManager")
#BiocManager::install("GEOquery")
#library("GEOquery")
# load data
#gse<- getGEO(file = "GSE46449_series_matrix.txt.gz", getGPL = FALSE) # miRNA expression profiling of human immune cell subsets (HUG)
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
# gse<- getGEO("GSE28489", destdir = ".", getGPL = FALSE)
gse <- gse[[1]]
# inspect matrix of expression values
ex <- exprs(gse)
dim(ex)
colnames(ex)






ex <- ex_TEST



# take log - not always necessary
ex2 <- log2(ex)
# analyze value distributions
boxplot(ex2)
ex2 <- na.omit(as.matrix(ex2))
# PCA
pca <- prcomp(t(ex2))
summary(pca)
screeplot(pca, main="PCA for ex_TEST")
# draw PCA plot
grpcol <- c(rep("blue",5), rep("red",5), rep("green",5), rep("yellow",5) , rep("pink",3), rep("purple",5), rep("orange",5) )
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2, ex_TEST", type="p", pch=10)#, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.5, labels = colnames(ex2))

#We can also do the PCA on certain special set : only the Frontal Cortex for example or only a certain treatment
# take log - not always necessary
ex2 <- log2(ex_Dr)
# analyze value distributions
boxplot(ex2)
ex2 <- na.omit(as.matrix(ex2))
# PCA
pca <- prcomp(t(ex2))
summary(pca)
screeplot(pca, main = "PCA for Dorsal Raphe")
# draw PCA plot
grpcol <- c(rep("blue",5), rep("red",5), rep("green",5), rep("yellow",5) , rep("pink",3), rep("purple",5), rep("orange",5) )
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2, Dorsal Raphe", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.5, labels = colnames())














############################## WEEK 4 K-Means and Hierachical Clustering#################################
#install.packages("BiocManager")
#BiocManager::install("GEOquery")
#library("GEOquery")
#BiocManager::install("useful")
#library("useful")
# miRNA expression profiling of human immune cell subsets
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
gse <- gse[[1]]
na.omit(gse)
gse<-new("ExpressionSet", exprs=as.matrix(gse))
ex <- exprs(gse)
dim(ex)






# K-Means


label=colnames(ex2)




ex2 <- log2(ex_Fc[,c(-9,-1)])

#boxplot(ex2)
ex2 <- na.omit(as.matrix(ex2))
pca <- prcomp(t(ex2))
grpcol <- c(rep("blue",5), rep("red",5), rep("green",5), rep("yellow",5), rep("pink",3), rep("purple",5), rep("orange",5) )
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10 )#, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.5, labels = colnames(ex2))
k <- 4
kmeans_result <- kmeans(t(ex2), k)
table(kmeans_result$cluster)
# Note that the cluster visualization is not trivial (space of 7815 dimensions!). We use the
# plot function in the 'useful' package that performs a dimensionality reduction using PCA
plot(kmeans_result, data=t(ex2))
plot(kmeans_result, data=t(ex2)) + geom_text(aes(label=colnames(ex2)),hjust=0,vjust=0)
#text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.5, labels = '')








#HIERARCHICAL CLUSTERING EXAMPLE
#install.packages("BiocManager")
#BiocManager::install("GEOquery")
#library("GEOquery")
#BiocManager::install("useful")
#library("useful")
# miRNA expression profiling of human immune cell subsets
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
gse <- gse[[1]]
ex <- exprs(gse)
dim(ex)
ex2 <- log2(ex)
boxplot(ex2)



ex2 <- na.omit(as.matrix(ex_TEST))
pca <- prcomp(t(ex2))
grpcol <- c(rep("blue",5), rep("red",5), rep("green",5), rep("yellow",5), rep("pink",3), rep("purple",5), rep("orange",5) )
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca$data), cex=0.5, labels = colnames(ex2))
dist_matrix <- dist(t(ex2))
hc_result <- hclust(dist_matrix, method = "centroid") # try with different methods of linkage
                                                 # i.e. complete, single, average
k <- 3
groups <- cutree(hc_result, k=k)    #try with differents types of cuts
table(groups)
plot(hc_result, hang <- -1, labels= colnames(ex2))
rect.hclust(hc_result, k = 3, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

#remove the outliers for overall data
ex_corr <- ex[,-c(1,19,23,24,41)]


#Remove the outliers for the ex_TEST
ex2 <- ex2[,-7]
ex2 <- ex2[,-1]
colnames(ex2)

#Try with dissimilarity measure = correlation matrix ! 


#Computation of the correlation matrix : 
install.packages("corrplot")
source("http://www.sthda.com/upload/rquery_cormat.r")
#mydata <- ex2
corr_matrix <- rquery.cormat(ex_corr, type="full")
corr_matrix <- corr_matrix$r
dissim1<-1-corr_matrix
dist1<-as.dist(dissim1)

hc_result <- hclust(dist1, method = "single") 
k <- 3
groups <- cutree(hc_result, k=k)    
table(groups)
plot(hc_result, main="Cluster Dendrogram all dataset, correlation matrix")
rect.hclust(hc_result, k = 3, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) # red boxes to show groups

str(corr_matrix)
head(corr_matrix$r)














############################## WEEK 5 : Random Forest #################################
#install.packages("BiocManager")
# load data
#BiocManager::install("ALL"); library("ALL"); data(ALL)
# keep only 30 arrays JUST for computational convenience
colnames(gse) <- c("Fc_MDMA/VLX_1","Fc_MDMA/VLX_2","Fc_MDMA/VLX_3","Fc_MDMA/SAL_1","Fc_MDMA/SAL_2",
                  "Fc_MDMA/SAL_3","Fc_MDMA/SAL_4","Fc_SAL/VLX_1","Fc_SAL/VLX_2","Fc_SAL/VLX_3",
                  "Fc_SAL/VLX_4","Fc_SAL/SAL_1","Fc_SAL/SAL_2","Fc_SAL/SAL_3","Fc_SAL/SAL_4",
                  "Hp_MDMA/VLX_1","Hp_MDMA/VLX_2","Hp_MDMA/VLX_3","Hp_MDMA/VLX_4","Hp_MDMA/SAL_1",
                  "Hp_MDMA/SAL_2","Hp_MDMA/SAL_3","Hp_MDMA/SAL_4","Hp_SAL/VLX_1","Hp_SAL/VLX_2",
                  "Hp_SAL/VLX_3","Hp_SAL/VLX_4","Hp_SAL/SAL_1","Hp_SAL/SAL_2","Hp_SAL/SAL_3",
                  "Hp_SAL/SAL_4","Dr_MDMA/VLX_1","Dr_MDMA/VLX_2","Dr_MDMA/VLX_3","Dr_MDMA/VLX_4",
                  "Dr_MDMA/SAL_1","Dr_MDMA/SAL_3","Dr_MDMA/SAL_4","Dr_SAL/VLX_1","Dr_SAL/VLX_2",
                  "Dr_SAL/VLX_3","Dr_SAL/VLX_4","Dr_SAL/SAL_1","Dr_SAL/SAL_2","Dr_SAL/SAL_3","Dr_SAL/SAL_4")
colnames(gse)
gse_TEST <- gse[,c(1,2,3,16,17,18,19,32,33,34,35,12,13,14,15,28,29,30,31,43,44,45,46)]
colnames(gse_TEST)
#Then we remove the two outliers
gse_TEST <- gse_TEST[,-7]
gse_TEST <- gse_TEST[,-1]
colnames(gse_TEST)

e.mat <- 2^(exprs(gse_TEST))


e.mat <- 2^(exprs(gse)[,c(1:46)])
# Genefilter package is very useful to preprocess data
# here we remove genes whose value is not > 100 in at least 20% of the samples
# filter genes on raw scale, then return to log scale;
#BiocManager::install("genefilter"); 


library("genefilter");
ffun <- filterfun(pOverA(0.20,100))
t.fil <- genefilter(e.mat,ffun)
small.eset <- log2(e.mat[t.fil,])
dim(small.eset) # 22466 genes, 23 arrays (11 positive and 12 negative)
group <- c(rep('MDMAVLX',9),rep('SALSAL',12)) # classification, in order
# Build RF
#BiocManager::install(" randomForest"); 
#library(randomForest)
set.seed(1234)
print(date())
rf <- randomForest(x=t(small.eset), y=as.factor(group), ntree=1000)
print(date()) # it takes about 20 seconds
# a trivial test
predict(rf, t(small.eset[, 1:5]))
print(date()) #so that we can compare with the previous to know the computation time 

rf
plot(rf, main='Random Forest for ex_TEST') 

##DRAWING A HEATMAP#

# Look at variable importance
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing="TRUE")
plot(c(1:nrow(small.eset)),imp.temp[t],log='x',cex.main=1.5,
     xlab='gene rank',ylab='variable importance',cex.lab=1.5,
     pch=16,main='ALL subset results, ex_TEST')
# Get subset of expression values for 25 most 'important' genes
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order
t <- is.element(rownames(small.eset),gn.25)
sig.eset <- small.eset[t,] # matrix of expression values, not necessarily in order
## Make a heatmap, with group differences obvious on plot
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
colnames(sig.eset) <- group # This will label the heatmap columns
csc <- rep(hmcol[50],23)
csc[group=='T'] <- hmcol[200]
# column side color will be purple for T and orange for B
heatmap(sig.eset, scale="row", col=hmcol, main ='Heatmap for ex_TEST dataset', ColSideColors=csc)













########### Random Forest : For the overall data ###########"
colnames(gse) <- c("Fc_MDMA/VLX_1","Fc_MDMA/VLX_2","Fc_MDMA/VLX_3","Fc_MDMA/SAL_1","Fc_MDMA/SAL_2",
                   "Fc_MDMA/SAL_3","Fc_MDMA/SAL_4","Fc_SAL/VLX_1","Fc_SAL/VLX_2","Fc_SAL/VLX_3",
                   "Fc_SAL/VLX_4","Fc_SAL/SAL_1","Fc_SAL/SAL_2","Fc_SAL/SAL_3","Fc_SAL/SAL_4",
                   "Hp_MDMA/VLX_1","Hp_MDMA/VLX_2","Hp_MDMA/VLX_3","Hp_MDMA/VLX_4","Hp_MDMA/SAL_1",
                   "Hp_MDMA/SAL_2","Hp_MDMA/SAL_3","Hp_MDMA/SAL_4","Hp_SAL/VLX_1","Hp_SAL/VLX_2",
                   "Hp_SAL/VLX_3","Hp_SAL/VLX_4","Hp_SAL/SAL_1","Hp_SAL/SAL_2","Hp_SAL/SAL_3",
                   "Hp_SAL/SAL_4","Dr_MDMA/VLX_1","Dr_MDMA/VLX_2","Dr_MDMA/VLX_3","Dr_MDMA/VLX_4",
                   "Dr_MDMA/SAL_1","Dr_MDMA/SAL_3","Dr_MDMA/SAL_4","Dr_SAL/VLX_1","Dr_SAL/VLX_2",
                   "Dr_SAL/VLX_3","Dr_SAL/VLX_4","Dr_SAL/SAL_1","Dr_SAL/SAL_2","Dr_SAL/SAL_3","Dr_SAL/SAL_4")
colnames(gse)
gse_TEST <- gse[,c(1,2,3,16,17,18,19,32,33,34,35,12,13,14,15,28,29,30,31,43,44,45,46)]
colnames(gse_TEST)
#Then we remove the two outliers
# gse_TEST <- gse_TEST[,-7]
# gse_TEST <- gse_TEST[,-1]
# colnames(gse_TEST)

# e.mat <- 2^(exprs(gse_TEST))


e.mat <- 2^(exprs(gse)[,c(1:46)])
# Genefilter package is very useful to preprocess data
# here we remove genes whose value is not > 100 in at least 20% of the samples
# filter genes on raw scale, then return to log scale;
#BiocManager::install("genefilter"); 


#library("genefilter");
ffun <- filterfun(pOverA(0.20,100))
t.fil <- genefilter(e.mat,ffun)
small.eset <- log2(e.mat[t.fil,])
dim(small.eset) # 22497 genes, 30 arrays (15 B and 15 T)
group <- c(rep('MDMAVLX',3),rep('MDMDASAL',4),rep('SALVLX',4),rep('SALSAL',4),
           rep('MDMAVLX',4),rep('MDMDASAL',4),rep('SALVLX',4),rep('SALSAL',4),
           rep('MDMAVLX',4),rep('MDMDASAL',3),rep('SALVLX',4),rep('SALSAL',4))# classification, in order

           
           
           # Build RF
#BiocManager::install(" randomForest"); 
#library(randomForest)
set.seed(1234)
print(date())
rf <- randomForest(x=t(small.eset), y=as.factor(group), ntree=1000)
print(date()) # it takes about 20 seconds
# a trivial test
predict(rf, t(small.eset))#[, 1:5]))
print(date()) #so that we can compare with the previous to know the computation time 

rf
plot(rf, main='Random Forest for all data')

##DRAWING A HEATMAP#

# Look at variable importance
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing="TRUE")
plot(c(1:nrow(small.eset)),imp.temp[t],log='x',cex.main=1.5,
     xlab='gene rank',ylab='variable importance',cex.lab=1.5,
     pch=16,main='ALL subset results')
# Get subset of expression values for 25 most 'important' genes
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order
t <- is.element(rownames(small.eset),gn.25)
sig.eset <- small.eset[t,] # matrix of expression values, not necessarily in order
## Make a heatmap, with group differences obvious on plot
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
colnames(sig.eset) <- group # This will label the heatmap columns
csc <- rep(hmcol[50],46)
csc[group=='T'] <- hmcol[200]
# column side color will be purple for T and orange for B
heatmap(sig.eset, scale="row", col=hmcol, main ='Heatmap for the overall data', ColSideColors=csc)



####### Random Forest : only for specific zones 

e.mat <- 2^(exprs(gse_TEST))
e.mat <- 2^(exprs(gse_Dr))

e.mat <- 2^(exprs(gse)[,c(1:46)])
# Genefilter package is very useful to preprocess data
# here we remove genes whose value is not > 100 in at least 20% of the samples
# filter genes on raw scale, then return to log scale;
#BiocManager::install("genefilter"); 


#library("genefilter");
ffun <- filterfun(pOverA(0.20,100))
t.fil <- genefilter(e.mat,ffun)
small.eset <- log2(e.mat[t.fil,])
dim(small.eset) # 22497 genes, 30 arrays (15 B and 15 T)
group <- c(#rep('MDMAVLX',3),rep('MDMDASAL',4),rep('SALVLX',4),rep('SALSAL',4)),
           #rep('MDMAVLX',4),rep('MDMDASAL',4),rep('SALVLX',4),rep('SALSAL',4),
           rep('MDMAVLX',4),rep('MDMDASAL',3),rep('SALVLX',4),rep('SALSAL',4))# classification, in order



# Build RF
#BiocManager::install(" randomForest"); 
#library(randomForest)
set.seed(1234)
print(date())
rf <- randomForest(x=t(small.eset), y=as.factor(group), ntree=1000)
print(date()) # it takes about 20 seconds
# a trivial test
predict(rf, t(small.eset))#[, 1:5]))
print(date()) #so that we can compare with the previous to know the computation time 

rf
plot(rf, main='Random Forest Dorsal Raphe')

##DRAWING A HEATMAP#

# Look at variable importance
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing="TRUE")
plot(c(1:nrow(small.eset)),imp.temp[t],log='x',cex.main=1.5,
     xlab='gene rank',ylab='variable importance',cex.lab=1.5,
     pch=16,main='ALL subset results, Frontal Cortex')
# Get subset of expression values for 25 most 'important' genes
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order
t <- is.element(rownames(small.eset),gn.25)
sig.eset <- small.eset[t,] # matrix of expression values, not necessarily in order
## Make a heatmap, with group differences obvious on plot
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
colnames(sig.eset) <- group # This will label the heatmap columns
csc <- rep(hmcol[50],15)
csc[group=='T'] <- hmcol[200]
# column side color will be purple for T and orange for B
heatmap(sig.eset, scale="row", col=hmcol, main ='Heatmap for Frontal Cortex', ColSideColors=csc)


############################## WEEK 6: LDA ############################
###### LDA 
install.packages("BiocManager")
BiocManager::install("genefilter"); library("genefilter")
BiocManager::install("GEOquery");
library("GEOquery")
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
length(gse)
gse <- gse[[1]]
show(gse)
ex <- exprs(gse)
dim(ex)
colnames(ex)







# ex2 <- log2(ex)
ex2<-log2(na.omit(ex)+40)
ex3<-ex2[,1:40]
library("genefilter")
f <- factor(c(rep(0,20), rep(1,20)))
tt40<-rowttests(ex3,f)
keepers<-which(p.adjust(tt40$p.value)<0.1)

tex3 <- t(ex3)
tex3 <- tex3[,keepers]
dat <- cbind(tex3,c(rep(0,20),rep(1,20)))
colnames(dat)[ncol(dat)] <- "AFFECTED"

n.controls<-10
n.affected<-10
train <- sample(1:(n.controls), (n.controls-5)) #we take an array of 20-5=15 samples for the training
test <- setdiff(1:(n.controls),train)
test<- c(test, test+20)
train <- c(train, train+20)

library("MASS")
mod <- lda(AFFECTED ~ ., data=as.data.frame(dat), prior =
             c(0.5,0.5), subset = train)
plot(mod)
mod.values <- predict(mod, as.data.frame(dat[train,]))
plot(mod.values$x[,1], ylab=c("LDA Axis"), type = "p", main = "LDA for Fc")
text(mod.values$x[,1],
     col=c(as.numeric(dat[train,"AFFECTED"])+10),
     labels = colnames(ex2),
     cex = 0.5)
preds<-predict(mod, as.data.frame(dat[test,]))
preds$class

########## ROC CURVE 
#install.packages("BiocManager")
#BiocManager::install("genefilter");
library("genefilter");
#BiocManager::install("GEOquery");
library("GEOquery")

library("MASS")
mod <- lda(AFFECTED ~ ., data=as.data.frame(dat), prior =
             c(0.5,0.5), subset = train)
plot(mod)
mod.values <- predict(mod, as.data.frame(dat[train,]))
mod.values$class
plot(mod.values$x[,1], ylab=c("LDA Axis"))
text(mod.values$x[,1],
     col=c(as.numeric(dat[train,"AFFECTED"])+10))
preds<-predict(mod, as.data.frame(dat[test,]))
preds$class
table(as.numeric(preds$class),
      
      as.numeric(dat[test,AFFECTED]) )

library("pROC")
roc_lda <- plot.roc(as.numeric(preds$class),
                    
                    as.numeric(dat[test,AFFECTED]) )

plot(roc_lda, col ="grey")


########## Version for ex_TEST for LDA 


# ex2 <- log2(ex)
ex2<-log2(na.omit(ex_TEST[,-15])+22)
ex3<-ex2[,1:22]
library("genefilter")
f <- factor(c(rep(1,11), rep(0,11)))
tt40<-rowttests(ex3,f)
keepers<-which((tt40$p.value)<0.05) #was at 0.1 initially

tex3 <- t(ex3)
tex3 <- tex3[,keepers]
dat <- cbind(tex3,c(rep(1,11), rep(0,11)))
colnames(dat)[ncol(dat)] <- "AFFECTED"
dim(dat)

n.controls<-11
n.affected<-11
train <- sample(1:(n.controls), (n.controls-5)) 
test <- setdiff(1:(n.controls),train)
test<- c(test, test+22)
train <- c(train, train+22)

library("MASS")
mod <- lda(AFFECTED ~ ., data=as.data.frame(dat), 
           prior =c(1/2,1/2), subset = train)
             
plot(mod)
mod.values <- predict(mod, as.data.frame(dat[train,]))
plot(mod.values$x[,1], ylab=c("LDA Axis"), type = "p", main = "LDA for ex_TEST")
text(mod.values$x[,1],
     col=c(as.numeric(dat[train,"AFFECTED"])+11),
     labels = colnames(ex2),
     cex = 0.5)
preds<-predict(mod, as.data.frame(dat[test,]))
preds$class

########## ROC CURVE 
#install.packages("BiocManager")
BiocManager::install("genefilter"); library("genefilter");
BiocManager::install("GEOquery");library("GEOquery")

library("MASS")
mod <- lda(AFFECTED ~ ., data=as.data.frame(dat), prior =
             c(0.5,0.5), subset = train)
plot(mod)
mod.values <- predict(mod, as.data.frame(dat[train,]))
mod.values$class
plot(mod.values$x[,1], ylab=c("LDA Axis"))
text(mod.values$x[,1],
     col=c(as.numeric(dat[train,"AFFECTED"])+10))
preds<-predict(mod, as.data.frame(dat[test,]))
preds$class
table(as.numeric(preds$class),
      
      as.numeric(dat[test,AFFECTED]) )

library("pROC")
roc_lda <- plot.roc(as.numeric(preds$class),
                    
                    as.numeric(dat[test,AFFECTED]) )

plot(roc_lda, col ="grey")

################### CARET ################## Used to simplify calculus by integrating pre-programmed functions
#source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
library("genefilter")
biocLite("GEOquery")
library("GEOquery")
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
length(gse)
gse <- gse[[1]]
show(gse)
ex <- exprs(gse)
dim(ex)
colnames(ex)


# ex2 <- log2(ex)
ex2<-log2(na.omit(ex_TEST[,-15])+22)
ex3<-ex2[,1:22]
library("genefilter")
f <- factor(c(rep(1,11), rep(0,11)))
tt40<-rowttests(ex3,f)
keepers<-which((tt40$p.value)<0.1)
tex3 <- t(ex3)
tex3 <- tex3[,keepers]
dat <- cbind(tex3,c(rep(1,11), rep(0,11)))
colnames(dat)[ncol(dat)] <- "AFFECTED"
#install.packages("caret")
library("caret")
inTrain <- createDataPartition(dat, p = .75, list = FALSE) #we take 75% of data for training and 25% for testing
train <- dat[ inTrain, ]
test <- dat[ -inTrain, ]
Ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

ldaFit <- train ( AFFECTED ~ .,
                  data = as.data.frame(train) ,
                  method = ' lda ' ,
                  trControl = ctrl ,
                  metric =  'ROC' )

plot(ldaFit, main ='LDA', labels=colnames(dat))
pred <- predict(ldaFit, newdata = test)

#We try to modify the very last paragraph : for the train we take the Fc and for the test the Hp
# ex2 <- log2(ex)
#ex2<-log2(na.omit(ex)+20)
#ex3<-ex2[,1:31]
#library("genefilter")




f <- factor(c(rep(0,10), rep(1,10)))
tt40<-rowttests(ex3,f)
keepers<-which(p.adjust(tt40$p.value)<0.1)
tex3 <- t(ex3) #take the transpose (interversion of lines and columns)
tex3 <- tex3[,keepers]
dat <- cbind(tex3,c(rep(0,10),rep(1,10)))
colnames(dat)[ncol(dat)] <- "AFFECTED"
#install.packages("caret")
library("caret")
inTrain <- createDataPartition(y = dat, p = .75, list = FALSE) #we take 75% of data for training and 25% for testing
train <- dat[ 1:inTrain, ] #problem of out of bounds
test <- dat[ inTrain+1:nrow(dat), ] # problem of out of bounds
Ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

ldaFit <- train ( AFFECTED ~ .,
                  data = as.data.frame(train) ,
                  method = ' lda ' ,
                  trControl = ctrl ,
                  metric =  'ROC' )

plot(ldaFit, main ='LDA', labels=colnames(dat))
pred <- predict(ldaFit, newdata = test)






#We do it by hand for the ex_TEST :  


# ex2 <- log2(ex)
#First we remove 3 samples out of the 23 in our data ex_TEST : 
ex_TEST2 <- ex_TEST[,-7]
colnames(ex_TEST2)
ex_TEST2 <- ex_TEST2[,-14]
colnames(ex_TEST2)
ex_TEST2 <- ex_TEST2[,-20]
colnames(ex_TEST2)

ex2<-log2(na.omit(ex_TEST2)+20)
ex3<-ex2 #[,1:20]
#library("genefilter")
f <- factor(c(rep(1,10), rep(0,10)))
tt40<-rowttests(ex3,f)
keepers<-which(p.adjust(tt40$p.value)<0.1) #select the  columns (supervised mode so p<0.1)
tex3 <- t(ex3) #take the transpose (interversion of lines and columns)
tex3 <- tex3[,keepers]
dat <- cbind(tex3,c(rep(1,10),rep(0,10))) #11 first are positive and the 12 others are negative
colnames(dat)[ncol(dat)] <- "AFFECTED"
#install.packages("caret")
#library("caret")
#inTrain <- createDataPartition(y = dat, p = .75, list = FALSE) #we take 75% of data for training and 25% for testing
train <- dat[ c(1:5,11:15), ] #problem of out of bounds
test <- dat[ c(6:10,16:20), ] # problem of out of bounds
Ctrl <- trainControl(method = "repeatedcv",
                     repeats = 3,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

ldaFit <- train ( AFFECTED ~ .,
                  data = as.data.frame(train) ,
                  method = ' lda ' ,
                  trControl = ctrl ,
                  metric =  'ROC' )

plot(ldaFit, main ='LDA', labels=colnames(dat))
pred <- predict(ldaFit, newdata = test)

####################### COMPARE LDA AND PCA ON THE SAME PLOT####################""""
require(MASS)
require(ggplot2)
require(scales)
require(gridExtra)


#We transpose and then add the AFFECTED column
# ex_TEST <- ex_TEST[,-7]
# ex_TEST <- ex_TEST[,-1]
# colnames(ex_TEST)
# tex_TEST <- t(ex_TEST)
# dim(tex_TEST)
# rownames(tex_TEST)


ex3 <- ex_TEST
f <- factor(c(rep(1,9), rep(0,12)))
tt40<-rowttests(ex3,f)
keepers<-which(p.adjust(tt40$p.value)<0.1) #select the  columns (supervised mode so p<0.1)
keepers<-which((tt40$p.value)<0.1) 
tex3 <- t(ex3) #take the transpose (interversion of lines and columns)
tex3 <- tex3[,keepers]

dat <- cbind(tex3,c(rep(1,9),rep(0,12))) #we put the value 0 times 10 and the value 1 times 10
colnames(dat)[ncol(dat)] <- "AFFECTED"
dim(dat)


pca <- prcomp(iris[,-5],
              center = TRUE,
              scale. = TRUE)

prop.pca = pca$sdev^2/sum(pca$sdev^2)

# lda <- lda(AFFECTED ~ ., 
#            as.data.frame(dat))

prop.lda = r$svd^2/sum(r$svd^2)

plda <- predict(object = lda,
                newdata = iris)

dataset = data.frame(species = iris[,"Species"],
                     pca = pca$x, lda = plda$x)

p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = species, shape = species), size = 2.5) + 
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))

p2 <- ggplot(dataset) + geom_point(aes(pca.PC1, pca.PC2, colour = species, shape = species), size = 2.5) +
  labs(x = paste("PC1 (", percent(prop.pca[1]), ")", sep=""),
       y = paste("PC2 (", percent(prop.pca[2]), ")", sep=""))

grid.arrange(p1, p2)















############################ WEEK 7 : Lasso  ###########################
##########LASSO
#install.packages("BiocManager")
BiocManager::install("GEOquery")
library("GEOquery")
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
length(gse)
# ... shows us there is only one object in it. We assign it to
# the same variable.
#gse <- gse[[1]]
# show what we have:
show(gse)
na.omit(gse)
# The actual expression data are accessible in the "exprs" ... an
# Expression Set, the generic data class that BioConductor uses
# for expression data
gse<-new("ExpressionSet", exprs=as.matrix(gse)) #replaces the exprs() function
head(exprs(gse))
length(gse)
gse <- gse[[1]]
show(gse)
ex <- exprs(gse)
dim(ex)
colnames(ex)











ex2<-log2(na.omit(ex)+40)
ex3<-ex2[,1:40]
dat <- t(ex3)
y <- c(rep(0,20),rep(1,20))
f <- factor(y)

install.packages("glmnet")
library("glmnet")
fit=glmnet(dat,y,standardize=FALSE,family="binomial")
plot(fit, xvar = "lambda", label=TRUE)
cfit=cv.glmnet(dat,y,standardize=FALSE,family="binomial")
plot(cfit)
coef(cfit, s=cfit$lambda.min)
# repeat analysis but by using train + test sample subsets
n.controls<-20
n.affected<-20
#Lines 339 to 347 : Non deterministic ! Try executing multiple times and compare the results 
train <- sample(1:(n.controls), (n.controls-5))
test <- setdiff(1:(n.controls),train)
test<- c(test, test+20)
train <- c(train, train+20)
fit=glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(fit)
cfit=cv.glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(cfit)
predict(fit,dat[test,], type="class", s= cfit$lambda.min)
# plot ROCR curve
library("ROCR")
pred2 <- predict(fit,dat[test,], type="response", s=cfit$lambda.min)
plot(performance(prediction(pred2, y[test]), 'tpr', 'fpr'))
# compute Area Under the Curve (AUC)
auc.tmp <- performance(prediction(pred2, y[test]),"auc")
auc <- as.numeric(auc.tmp@y.values)
print(auc)

########### CARET for Lasso
#install.packages("BiocManager")
BiocManager::install("GEOquery")
library("GEOquery")
gse<- getGEO(file = "GSE47541_series_matrix.txt", getGPL = FALSE)
length(gse)
# ... shows us there is only one object in it. We assign it to
# the same variable.
#gse <- gse[[1]]
# show what we have:
show(gse)
na.omit(gse)
# The actual expression data are accessible in the "exprs" ... an
# Expression Set, the generic data class that BioConductor uses
# for expression data
gse<-new("ExpressionSet", exprs=as.matrix(gse)) #replaces the exprs() function
head(exprs(gse))
length(gse)
gse <- gse[[1]]
show(gse)
ex <- exprs(gse)
dim(ex)
colnames(ex)









ex2<-log2(na.omit(ex)+30)
ex3<-ex2[,1:30]
y <- c(rep('MDMAVLX',3),rep('MDMDASAL',4),rep('SALVLX',4),rep('SALSAL',4),
        rep('MDMAVLX',4),rep('MDMDASAL',4),rep('SALVLX',4),rep('SALSAL',3))
# IMPORTANT: to describe outcome you need to:
# - use factor (will force classification model)
# - with labels (caret wants proper class IDs)
#f <- factor(y, labels = c("Control", "Affected"))
dat <- t(ex3)

#install.packages("pROC")
library("pROC")
#install.packages("caret")
library("caret")
inTrain <- createDataPartition(y, p = .75, list = FALSE)
train <- dat[ inTrain, ]
test <- dat[ -inTrain, ]
ctrl <- trainControl(method = "cv",
                     
                     classProbs = TRUE,
                     summaryFunction = fourClassSummary)

lassoFit <- train ( train,
                    f[inTrain],
                    
                    method = "glmnet",
                    family = "binomial",
                    tuneGrid = expand.grid(alpha = 1,
                                           lambda = seq(0,1,by=0.05)),
                    
                    trControl = ctrl,
                    metric = "ROC" )

plot(lassoFit)
pred <- predict(lassoFit, newdata = test, type="prob")
pred2 <- predict(lassoFit, test, type="raw",
                 s=lassoFit$finalModel$lambdaOpt)

auc <- roc(f[-inTrain], pred[[2]])
print(auc$auc)
plot(auc)

######### Week 8 = SCUDO ####################
######### Weeks 9 & 10 = DAVID ##############
#Conversion of the genes_ID
# table = read.csv(file = "GEO_Accession_viewer_GSE47541.txt", sep = "\t")
# # Create a small program that takes as input the old list of imporant genes and, knowing the TABLE 
# # then we can do the correspandance and have access to the good expression of genes
# # Or either use the lumirat package directly to get the correspondance
# #Install the package, then get the lumiratAllSymbol and execute the first 7 lines of codes in the descirption
BiocManager::install("lumiRatAll.db")
library("lumiRatAll.db")
x <- lumiRatAllSYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five probes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

###### Correspondance when table 1 contains the signature's ID in the constructor syntax ! 

table1 <- read.csv(file="SCUDO_testing_signatures(overalldata-ColumnsNotWellAssigned).txt", sep = "\t") # liste des signatures 
table2 <- descrip #C'est le dictionnaire contenant les colonnes ID en 1 puis ILM_XXX en 2 
colnames(table1)
colnames(table2)

#get the signatures : do the correspondance between these columns, we'd like to have a list with
# a list of the ILM_XXXXX signatures genes as output
head(table1[[17]]) #signatures i constructor's ID 
head(table2[[2]]) #ILM_XXX
head(table2[[1]]) #Constructor ID's
head(table2[[3]]) #Symbol



#install.packages("dplyr")
#require(dplyr)
  
table1 <- read.table(file ="SCUDO_testing_signatures_MDMA_Dr.sig.txt", skip = 2, nrows = 250)


DF=data.frame(X=c(table1[9])) # list of constructor's ID that we want to convert 
Dico=data.frame(X=c(table2[1]),Y=c(table2[2])) # il faut que Dico soit au format "data.frame" et non "list"

SIGN_ILM <- Dico$Target..[match(DF[,1], Dico$X.ID.., nomatch = NA)]
SIGN_ILM <- na.omit(SIGN_ILM)
SIGN_ILM <- as.data.frame(SIGN_ILM)
head(SIGN_ILM)
dim(SIGN_ILM)
SIGN_ILM <- SIGN_ILM[c(-1,-2),]
write.table(x = SIGN_ILM, file = "SIGN_ILM_MDMA_Dr_UP.txt", row.names = FALSE, col.names = FALSE, na = '', quote = FALSE)

##### Dictionnaire qui retourne les Symbol et prend en entrée les ID constructor

table1 <- read.table(file ="SCUDO_testing_signatures_VLX_Dr-DOWN.sig.txt", skip = 2, nrows = 500)
dim(table1)

DF=data.frame(X=c(table1[10]))
head(DF)
dim(DF)# list of constructor's ID that we want to convert 
Dico=data.frame(X=c(table2[1]),Y=c(table2[3])) # il faut que Dico soit au format "data.frame" et non "list"

SIGN_ILM <- Dico$Symbol..[match(DF[,1], Dico$X.ID.., nomatch = NA)]
SIGN_ILM <- na.omit(SIGN_ILM)
SIGN_ILM <- as.data.frame(SIGN_ILM)
head(SIGN_ILM)
dim(SIGN_ILM)
SIGN_ILM <- SIGN_ILM[c(-1,-2),]
# enlever les LOC 
# enelever juste le suffixe 
write.table(x = SIGN_ILM, file = "SIGN_Symbol_VLX_Dr_DOWN.txt", row.names = FALSE, col.names = FALSE, na = '', quote = FALSE)

####### Dictionnaire qui retourne les Symbol et prend en entrée les ILM_XXX
signatures <- read.table(file = "SIGN_ILM_MDMASAL.txt")
dim(signatures)
head(signatures)

DF=signatures # list of constructor's ID that we want to convert 
Dico=data.frame(X=c(table2[3]),Y=c(table2[2])) # il faut que Dico soit au format "data.frame" et non "list"

SIGN_ILM <- Dico$Symbol..[match(DF[,1], Dico$Target.., nomatch = NA)]
SIGN_ILM <- na.omit(SIGN_ILM)
SIGN_ILM <- as.data.frame(SIGN_ILM)
head(SIGN_ILM)
dim(SIGN_ILM)
SIGN_ILM <- SIGN_ILM[c(-1,-2),]
write.table(x = SIGN_ILM, file = "SIGN_Symbol_MDMASAL.txt", row.names = FALSE, col.names = FALSE, na = '', quote = FALSE)



#########????????? WEEK 11&12 = Cytoscape #######
#-> Do the SCUDO analysis for ex_MDMA, ex_VLX, ex_TEST then input into Cytoscape 