args <- commandArgs(T)


Data.File <- trimws(args[1])
Pheno.File <- trimws(args[2])
Covars.Factor <- trimws(args[3])
Covars.Num <- trimws(args[4])
OutPrefix <- trimws(args[5])

cat("Input arguments:\n")
cat("   Data file: ",Data.File , "\n")
cat("   Phenotype information file: ",Pheno.File , "\n")
cat("   Factor Covariates: ",Covars.Factor , "\n")
cat("   Numeric Covariates: ",Covars.Num , "\n")
cat("   Output files prefix: ",OutPrefix , "\n")

cat("\n")
cat("Loading libraries...\n")

suppressMessages(library(stringr))
suppressMessages(library(WGCNA))
suppressMessages(library(ggplot2))

cat("\n")

cat("Reading the data...\n")
Data = readRDS(Data.File)
Pheno = read.csv(Pheno.File,stringsAsFactors = F,row.names = 1)
Covars.Num = str_split(Covars.Num,pattern = ',',simplify = T)[1,]
Covars.Factor = str_split(Covars.Factor,pattern = ',',simplify = T)[1,]

for (i in 1:length(Covars.Factor)) {
  Covars.Factor[i] <- trimws(Covars.Factor[i])
  Pheno[,Covars.Factor[i]] <- as.factor(Pheno[,Covars.Factor[i]])
}
for (i in 1:length(Covars.Num)) {
  Covars.Num[i] <- trimws(Covars.Num[i])
  Pheno[,Covars.Num[i]] <- as.numeric(Pheno[,Covars.Num[i]])
}

all_var = c(Covars.Factor,Covars.Num)

cat("Calculate PCs...\n")
pc.object <- prcomp(t(Data))
pc <- pc.object$x
pc.importance <- summary(pc.object)$importance[2,]

cor_ <- matrix(data = NA,nrow =10,ncol = length(all_var) )
colnames(cor_) <- c(Covars.Factor , Covars.Num)
rownames(cor_) <- paste0(colnames(pc)[1:10]," (",round(pc.importance[1:10],digits = 2),"%)")
cor_pval <- cor_

cat("\n")

cat("Calculatin correlations...\n")
for (i in 1:10) {
  for (j in 1:length(all_var)) {
    var_ <- all_var[j]
    
    if(var_ %in% Covars.Factor){
      res1<-cor.test(as.numeric(pc[,i]),as.numeric(as.factor(Pheno[,var_])), method="spearman",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }else{
      res1<-cor.test(as.numeric(pc[,i]),as.numeric(Pheno[,var_]), method="pearson",exact = FALSE)
      cor_[i,var_]<-as.numeric(res1$estimate)
      cor_pval[i,var_]<-as.numeric(res1$p.value)
    }
  }
}

textMatrix = paste(signif(cor_, 2), "\n(",signif(cor_pval, 1), ")", sep = "")
dim(textMatrix) = dim(cor_)

cat("\n")

cat("Saving Plots...\n")

pdf(file =paste0( OutPrefix,".BatchDetection.pdf"))
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = cor_,
               xLabels = colnames(cor_),
               yLabels = rownames(cor_),
               ySymbols = rownames(cor_),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("PCA Analysis"))


pc=cbind.data.frame(pc,Pheno)

for (i in 1:length(Covars.Num)) {
  p1 <- ggplot(data = pc,aes_string(x = "PC1",y = "PC2",colour=Covars.Num[i]))+
    geom_point()+ xlab(rownames(cor_)[1]) + ylab(rownames(cor_)[2]) +ggtitle(Covars.Num[i])+
    scale_color_gradientn(colours = rainbow(5)) 
  print(p1)
}

Data.ColMeans <- colMeans(Data)
Data.ColMeans <- cbind.data.frame(Data.ColMeans , Pheno[,Covars.Factor])

for (i in 1:length(Covars.Factor)) {
  p1 <- ggplot(data = pc,aes_string(x = "PC1",y = "PC2",colour=Covars.Factor[i]))+
    geom_point() + xlab(rownames(cor_)[1]) + ylab(rownames(cor_)[2])+
    ggtitle(Covars.Factor[i]) 
  
  p2 <- ggplot() + 
    geom_boxplot(data = Data.ColMeans, aes_string(x=Covars.Factor[i],
                                                  y="Data.ColMeans",fill=Covars.Factor[i])) +
    ggtitle(Covars.Factor[i]) + ylab("Data Average Value for each Sample")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  p3 <- ggplot() + 
    geom_boxplot(data = pc, aes_string(x=Covars.Factor[i],
                                                  y="PC1",fill=Covars.Factor[i])) +
    ggtitle(Covars.Factor[i]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  p4 <- ggplot() + 
    geom_boxplot(data = pc, aes_string(x=Covars.Factor[i],
                                       y="PC2",fill=Covars.Factor[i])) +
    ggtitle(Covars.Factor[i]) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  
  if(length(levels(Data.ColMeans[,Covars.Factor[i]])) > 5){
    p1 <- p1 + theme(legend.position="none")
    p2 <- p2 + theme(legend.position="none")
    p3 <- p3 + theme(legend.position="none")
    p4 <- p4 + theme(legend.position="none")
  }
  print(p1)
  print(p2)
  print(p3)
  print(p4)
}

graphics.off()
