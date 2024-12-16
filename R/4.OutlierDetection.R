args <- commandArgs(TRUE)

expr_file <- trimws(args[1])
out_pref <- trimws(args[3])

cat("Input arguments:\n")
cat("    Expression data file:",expr_file,"\n")
cat("    Output prefix:",out_pref,"\n")
cat("\n")


cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dendextend))
suppressMessages(library(factoextra))
suppressMessages(library(Rtsne))

source("./mahalanobis.outlier.R")

cat("Reading the data...\n")
expr_ <- readRDS(file=expr_file)

cat("Finding outliers using Mahalanobis distance and Chi-squared dstribution...\n")
mahala.tsne <- mahalanobis.outlier(Data = expr_ , method = "tsne", seed=12345)
mahala.pca <- mahalanobis.outlier(Data = expr_ , method = "pca")

tiff(filename = paste0(out_pref,".Outlier.tSNE.Mahalanobis.tif"), res = 300, width = 6 , height = 6 , units = "in")
print(mahala.tsne$Plot.2D)
graphics.off()

tiff(filename = paste0(out_pref,".Outlier.PCA.Mahalanobis.tif"), res = 300, width = 6 , height = 6 , units = "in")
print(mahala.pca$Plot.2D)
graphics.off()

outliers <- data.frame(Sample = rownames(mahala.pca$Data.2D) , Outlier.PCA = mahala.pca$Data.2D$Outlier , 
                       Outlier.tSNE = mahala.tsne$Data.2D$Outlier)
write.csv(outliers , file = paste0(out_pref,"Outliers.csv"))
cat("hierarchical clustering...\n")

cat("  Euclidean distance...\n")
distance <- dist(t(expr_),method = "euclidean")
hclust.euc.av <- hclust(distance,method = "average")
hclust.euc.ward <- hclust(distance,method = "ward.D2")
hclust.euc.sing <- hclust(distance,method = "single")    #This method select min distance 
hclust.euc.comp <- hclust(distance,method = "complete")  #This method select max distance
hclust.euc.cent <- hclust(distance,method = "centroid")

cat("  Correlation-based distance...\n")
cols.cor <- cor(expr_, use = "pairwise.complete.obs", method = "pearson")
distance <- as.dist(1 - cols.cor)
hclust.cor.av <- hclust(distance,method = "average")
hclust.cor.ward <- hclust(distance,method = "ward.D2")
hclust.cor.sing <- hclust(distance,method = "single")
hclust.cor.comp <- hclust(distance,method = "complete")
hclust.cor.cent <- hclust(distance,method = "centroid")  

labelCol <- rep("black" , times = ncol(expr_))
labelCol[outliers$Outlier.PCA == "Yes"] = "blue"
labelCol[outliers$Outlier.tSNE == "Yes"] = "green"
labelCol[(outliers$Outlier.PCA == "Yes") & (outliers$Outlier.tSNE == "Yes")] = "red"

label.size = 0.7
title.size = 1.5
hclust.euc.av <-  hclust.euc.av %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size) 
hclust.euc.comp <-  hclust.euc.comp %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size) 
hclust.euc.ward <-  hclust.euc.ward %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size) 
hclust.euc.sing <-  hclust.euc.sing %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size)
hclust.euc.cent <-  hclust.euc.cent %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size)


hclust.cor.av <-  hclust.cor.av %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size) 
hclust.cor.comp <-  hclust.cor.comp %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size) 
hclust.cor.ward <-  hclust.cor.ward %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size) 
hclust.cor.sing <-  hclust.cor.sing %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size) 
hclust.cor.cent <-  hclust.cor.cent %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol) %>% set("leaves_col",labelCol) %>% 
  set("labels_cex", label.size)

tiff(filename = paste0(out_pref,"Outlier.hclust.tif"),res = 600 , units = "in" , width = 8 , height = 14)
par(mfrow=c(10,1), mar=c(2,3,2,1))
plot(hclust.euc.av, main = "Euclidean distance, Average method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.euc.sing, main = "Euclidean distance, Single method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.euc.comp, main = "Euclidean distance, Complete method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.euc.ward, main = "Euclidean distance, Ward.D2 method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.euc.cent, main = "Euclidean distance, Ward.D2 method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.cor.av, main = "Correlation distance, Average method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.cor.sing, main = "Correlation distance, Single method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.cor.comp, main = "Correlation distance, Complete method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.cor.ward, main = "Correlation distance, Ward.D2 method", xlab = "", sub = "", cex.main = title.size)
plot(hclust.cor.cent, main = "Correlation distance, Ward.D2 method", xlab = "", sub = "", cex.main = title.size)
graphics.off()
