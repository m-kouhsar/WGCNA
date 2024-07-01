args <- commandArgs(TRUE)

beta_file <- args[1] 
out_pref <- args[2]

cat("Input arguments:\n")
cat("    Beta value file:",beta_file,"\n")
cat("    Output prefix:",out_pref,"\n")
cat("\n")


cat("Loading libraries...\n")
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dendextend))

cat("Reading beta values...\n")
betas <- readRDS(file=beta_file)

## finding outlier samples
cat("Hclust...\n")

cat("  Euclidean distance...\n")
distance <- dist(t(betas),method = "euclidean")
hclust.euc.ave <- hclust(distance,method = "average")
hclust.euc.ward <- hclust(distance,method = "ward.D2")
hclust.euc.sing <- hclust(distance,method = "single")    #This method select min distance 
hclust.euc.comp <- hclust(distance,method = "complete")  #This method select max distance
hclust.euc.cent <- hclust(distance,method = "centroid")

cat("  Correlation-based distance...\n")
cols.cor <- cor(betas, use = "pairwise.complete.obs", method = "pearson")
distance <- as.dist(1 - cols.cor)
hclust.cor.ave <- hclust(distance,method = "average")
hclust.cor.ward <- hclust(distance,method = "ward.D2")
hclust.cor.sing <- hclust(distance,method = "single")
hclust.cor.comp <- hclust(distance,method = "complete")
hclust.cor.cent <- hclust(distance,method = "centroid")

cat("PCA...\n")
pc <- prcomp(t(betas))
pc.importance <- summary(pc)$importance[2,]
pc <- as.data.frame(pc$x[,1:2])

p1 <- ggplot() +geom_point(data = pc,aes(x = PC1,y = PC2))+
      xlab(paste0("PC1 ","(",round(pc.importance[1],digits = 2),")")) + 
      ylab(paste0("PC2 ","(",round(pc.importance[2],digits = 2),")")) + ggtitle("Samples PCA plot")

cat("Saving plots...\n")
pdf(file = paste0(out_pref,".OutlierDetection.hclust.pdf"),width = 24,height = 20)
par(cex = 0.6)
par(mar = c(0,4,2,0)) 

plot(hclust.cor.cent, main = "Sample clustering, distance=(1-correlation), method=Centroid", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.cor.comp, main = "Sample clustering, distance=(1-correlation), method=Complete", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.cor.sing, main = "Sample clustering, distance=(1-correlation), method=Single", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.cor.ward, main = "Sample clustering, distance=(1-correlation), method=Ward.D2", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.cor.ave, main = "Sample clustering, distance=(1-correlation), method=Average", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

plot(hclust.euc.cent, main = "Sample clustering, distance=euclidean, method=Centroid", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.euc.comp, main = "Sample clustering, distance=euclidean, method=Complete", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.euc.sing, main = "Sample clustering, distance=euclidean, method=Single", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.euc.ward, main = "Sample clustering, distance=euclidean, method=Ward.D2", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
plot(hclust.euc.ave, main = "Sample clustering, distance=euclidean, method=Average", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

graphics.off()

pdf(file = paste0(out_pref,".OutlierDetection.PCA.pdf"))
print(p1)

graphics.off()

cat("All done!")
cat("\n")


