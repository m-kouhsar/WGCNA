
args = commandArgs(T)

script_dir <- trimws(args[1])
expr.file <- trimws(args[2])
pheno.file <- trimws(args[3])
var.fact <- trimws(args[4])
var.num <- trimws(args[5])
OutPrefix <- trimws(args[6])

message("Input arguments:")
message("        Methylatio/Expression data file: ", expr.file)
message("        Phenotype file: ", pheno.file)
message("        Numeric variables: ", var.num)
message("        Categorical variables: ", var.fact)
message("        Output files prefix: ", OutPrefix)
cat("\n")

message("loading requiered packages...")
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(edgeR)))
suppressWarnings(suppressMessages(library(WGCNA)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(ggcorrplot)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(dendextend)))
cat("\n")
########################################################################
#
#          Reading the data
#
########################################################################
dir.create(dirname(OutPrefix) , recursive = T , showWarnings = F)
source(paste0(script_dir , "/R/DataChecking.Functions.R"))
message("Reading the data ...")
if(str_ends(string = expr.file , pattern = ".rds")){
  expr_mat <- readRDS(expr.file)
}else{
  expr_mat <- read.table(file = expr.file, stringsAsFactors = F,header = T, row.names = 1,check.names=F)
}

pheno <- read.csv(pheno.file , row.names = 1 , stringsAsFactors = F)

if(!identical(colnames(expr_mat) , rownames(pheno))){
  message("Warning message:\nColnames in the count matrix are not equal to the rownames in the phenotype file!\nShared names will be considered.")
  index <- intersect(colnames(expr_mat), rownames(pheno))
  message("Number of shared names: " , length(index))
  if(length(index) < 1){
    stop("There is no shared ID between Phneotype and Count data!")
  }
  expr_mat <- expr_mat[,index]
  pheno <- pheno[index , ]
}

var.fact <- trimws(str_split_1(var.fact , pattern = ","))
var.num <- trimws(str_split_1(var.num , pattern = ","))

pheno <- pheno[,c(var.fact,var.num)]

for (i in 1:length(var.fact)) {
  pheno[,var.fact[i]] <- as.factor(pheno[,var.fact[i]])
}
for (i in 1:length(var.num)) {
  pheno[,var.num[i]] <- as.numeric(pheno[,var.num[i]])
}

########################################################################
#
#          Checking for possible outliers and batch effects
#
########################################################################
plot_data <- cbind.data.frame(sample = colnames(expr_mat) ,pheno)
cor_plot1 <- generate_corr_analysis(data = subset(plot_data , select =c(-sample)),
                                   categorical_vars =var.fact ,
                                   plot.title = "Correlation test on all covariates" )
message("PCA analysis...")
pcs_obj <- prcomp(t(expr_mat), rank. = 10)
PCs <- as.data.frame(pcs_obj$x)

plot_data <- cbind.data.frame(plot_data , PCs)

plot_data$PC1_z <- scale(PCs$PC1)
plot_data$PC2_z <- scale(PCs$PC2)
plot_data$Outliers.PC.ZScore <- (abs(plot_data$PC1_z) > 3 | abs(plot_data$PC2_z) > 3)

plot_data$AveExpr <- colMeans(expr_mat , na.rm = T)
PC.PVar <- round(summary(pcs_obj)$importance[2,1:10]*100,digits = 2)

all_var = c(var.fact , var.num)

cor_plot2 <- generate_corr_analysis(data = subset(plot_data , select =c(-sample,-Outliers.PC.ZScore,-PC1_z , -PC2_z)),
                                   categorical_vars =var.fact ,
                                   plot.title = "Correlation test on all covariates" )

cor_ <- cor_plot2$cor.mat[paste0("PC",c(1:10)),all_var]
cor_pval <- cor_plot2$p.mat[paste0("PC",c(1:10)),all_var]
impact_val <- matrix(data = NA,nrow =10,ncol = length(all_var) )
colnames(impact_val) <- all_var
rownames(impact_val) <- colnames(PCs)[1:10]

for (i in 1:10) {
  for (j in 1:length(all_var)) {
    var_ <- all_var[j]
    
    if(cor_pval[i,var_] < 0.05){
      impact_val[i , var_] <- cor_[i,var_]**2 * PC.PVar[i]
    }else{
      impact_val[i , var_] = 0
    }
  }
}

textMatrix_cor = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix_cor) = dim(cor_)

pdf(file = paste0(OutPrefix , ".PCA.pdf"),width = 10,height = 10)

print(cor_plot1$plot)

par(mar = c(12, 8,3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = paste0(rownames(cor_),"(",PC.PVar[1:10],"%)"),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_cor,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "PCA Analysis (Correlation and P-value)")

textMatrix_imp = signif(impact_val, 2)
dim(textMatrix_imp) = dim(cor_)

par(mar = c(12, 8,3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = paste0(rownames(cor_),"(",PC.PVar[1:10],"%)"),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_imp,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               legendLabel = "Cor",
               main = expression("Significant PC Correlations" ~ (Impact ~ value == r^2 %*% "%Var"[PC])))


PCA.mahal <- mahalanobis.outlier(Data = expr_mat , method = "pca" , plot.title = "Outliers by Mahalanobis Distance")
plot_data$mdist <- PCA.mahal$Data.2D$mdist
plot_data$pchisq <- PCA.mahal$Data.2D$pchisq
plot_data$qchisq <- PCA.mahal$Data.2D$qchisq
plot_data$Outliers.Mahalanobis <- PCA.mahal$Data.2D$Outlier

print(PCA.mahal$Plot.2D)
print(PCA.mahal$Plot.QQ)

for (var_ in c(var.fact , var.num)) {
  p <- ggplot(data = plot_data, aes_string(x = "PC1_z" , y = "PC2_z" , colour = var_))  + 
    geom_point() + 
    # Add the +/- 3 z-score lines
    geom_vline(xintercept = 3, linetype = "dashed", color = "red") +
    geom_vline(xintercept = -3, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -3, linetype = "dashed", color = "red") +
    ggtitle("PCA Plot (Z-Scores) with Outlier Thresholds")+
    xlab("PC1 (Z Score)")+
    ylab("PC2 (Z Score)")+
    theme_bw()
  print(p)
  if(var_ %in% var.fact){
    p <- ggplot() + geom_boxplot(data = plot_data , aes_string(x = var_ ,  y = "AveExpr" , fill = var_))+
      ylab("Average gene expression")
    print(p)
  }
}

graphics.off()
######################################################################
message("Hirarchical clustering...")

cat("  Euclidean distance...\n")
distance <- dist(t(expr_mat) , method = "euclidean")
hclust.euc.av <- hclust(distance,method = "average")
hclust.euc.ward <- hclust(distance,method = "ward.D2")
hclust.euc.sing <- hclust(distance,method = "single")    #This method select min distance 
hclust.euc.comp <- hclust(distance,method = "complete")  #This method select max distance
hclust.euc.cent <- hclust(distance,method = "centroid")

cat("  Correlation-based distance...\n")
cols.cor <- cor(expr_mat, use = "pairwise.complete.obs", method = "pearson")
distance <- as.dist(1 - cols.cor)
hclust.cor.av <- hclust(distance,method = "average")
hclust.cor.ward <- hclust(distance,method = "ward.D2")
hclust.cor.sing <- hclust(distance,method = "single")
hclust.cor.comp <- hclust(distance,method = "complete")
hclust.cor.cent <- hclust(distance,method = "centroid")  

labelCol <- rep("black" , times = ncol(expr_mat))
labelCol[plot_data$Outliers.Mahalanobis == "Yes"] = "orange"
labelCol[plot_data$Outliers.PC.ZScore] = "coral"
labelCol[plot_data$Outliers.PC.ZScore & (plot_data$Outliers.Mahalanobis == "Yes")] = "red"

label.size = 0.7
title.size = 1.5

hclust.euc.av <-  hclust.euc.av %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size) 

hclust.euc.comp <-  hclust.euc.comp %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size) 

hclust.euc.ward <-  hclust.euc.ward %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size) 

hclust.euc.sing <-  hclust.euc.sing %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size)

hclust.euc.cent <-  hclust.euc.cent %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size)



hclust.cor.av <-  hclust.cor.av %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size) 

hclust.cor.comp <-  hclust.cor.comp %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size) 

hclust.cor.ward <-  hclust.cor.ward %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size) 

hclust.cor.sing <-  hclust.cor.sing %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size) 

hclust.cor.cent <-  hclust.cor.cent %>% as.dendrogram %>% hang.dendrogram %>%
  set("labels_colors",labelCol, order_value = TRUE) %>% 
  set("leaves_col",labelCol, order_value = TRUE) %>% 
  set("labels_cex", label.size)

pdf(file = paste0(OutPrefix , ".hClust.pdf"),width = 18,height = 30)
par(mfrow=c(10,1), mar=c(2,3,2,1))
plot(hclust.euc.av, main = "Euclidean distance, Average method", xlab = "", sub = "", cex.main = title.size)
legend("topright",
       legend = c("Mahalanobis", "PC Z-score", "Both"),
       col = c("orange", "pink", "red"),
       title = "Otliers",pch = c(16, 16))
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

#################################################
write.csv(plot_data , file = paste0(OutPrefix , ".DataChecking.csv") , row.names = F)

message("All done!")


