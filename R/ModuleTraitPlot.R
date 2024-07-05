data.frame.relationship <- function(df1 , df2, method="prsp", Plot=T,plot.title="",return_melt=T, categorical_columns){
  # method:
  #   pr -> Pearson correlation
  #   sp -> Spearman correlation
  #   lm -> Linear regression
  #   prsp -> Pearson correlation for numeric variables and Spearman correlation for categorical variables
  #   tsppr -> Pearson correlation for two numeric variables,
  #            Spearman correlation for two categorical variables or one categorical (more than two groups) and one numeric variables,
  #            t-test for a binary and a numeric variable
  
  method <- match.arg(method,c("pr","sp","prsp","lm","anovapr","tsppr"),several.ok = F)
  if(!is.data.frame(df1)){
    stop("df1 must be a dataframe")
  }
  if(!is.data.frame(df1)){
    stop("df2 must be a dataframe")
  }
  if(nrow(df1) != nrow(df2)){
    stop("Two dataframes must have the same number of rows")
  }
  cor_val <- matrix(data = NA, nrow = ncol(df1), ncol = ncol(df2))
  colnames(cor_val) <- names(df2)
  rownames(cor_val) <- names(df1)
  cor_p <- cor_val
  cor_plot <- NA
  
  if(method=="pr"){
    for(i in names(df1)){
      for (j in names(df2)) {
        res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="pearson",exact = FALSE)
        cor_val[i,j]<-res1$estimate
        cor_p[i,j]<-res1$p.value
      }
    }
  }
  if(method=="sp"){
    for(i in names(df1)){
      for (j in names(df2)) {
        res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="spearman",exact = FALSE)
        cor_val[i,j]<-res1$estimate
        cor_p[i,j]<-res1$p.value
      }
    }
  }
  if(method=="prsp"){
    for(i in names(df1)){
      for (j in names(df2)) {
        if((i %in% categorical_columns)|(j %in% categorical_columns)){
          res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="spearman",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }else{
          res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="pearson",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }
      }
    }
  }
  
  # https://stackoverflow.com/questions/52811684/running-a-two-sample-t-test-with-unequal-sample-size-in-r
  if(method=="tsppr"){
    for(i in names(df1)){
      for (j in names(df2)) {
        if(i %in% categorical_columns){
          flag1=T
        }else{
          flag1=F
        }
        if(j %in% categorical_columns){
          flag2=T
        }else{
          flag2=F
        }
        if(!flag1 & !flag2){
          res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="pearson",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }
        if(flag1 & !flag2){
          df <- cbind.data.frame(c1=df1[,i],c2=df2[,j])
          df.split <- split(df,as.factor(df$c1),drop = T)
          if(length(df.split) == 2){
            if(length(df.split[[1]]$c2) == length(df.split[[2]]$c2)){
              res1 <- t.test(df.split[[1]]$c2, df.split[[2]]$c2, paired = T)
              cor_val[i,j]<-res1$estimate
              cor_p[i,j]<-res1$p.value
            }else{
              if(var.test(df.split[[1]]$c2, df.split[[2]]$c2)$p.value < 0.05){
                res1 <- t.test(df.split[[1]]$c2, df.split[[2]]$c2, var.equal = F)
                cor_val[i,j]<-abs(res1$estimate[1] - res1$estimate[2])
                cor_p[i,j]<-res1$p.value
              }
              else{
                res1 <- t.test(df.split[[1]]$c2, df.split[[2]]$c2, var.equal = T)
                cor_val[i,j]<-abs(res1$estimate[1] - res1$estimate[2])
                cor_p[i,j]<-res1$p.value
              }
            }
          }else{
            res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="spearman",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }

        }
        if(!flag1 & flag2){
          df <- cbind.data.frame(c1=df1[,i],c2=df2[,j])
          df.split <- split(df,as.factor(df$c2),drop = T)
          if(length(df.split) == 2){
            if(length(df.split[[1]]$c1) == length(df.split[[2]]$c1)){
              res1 <- t.test(df.split[[1]]$c1, df.split[[2]]$c1, paired = T)
              cor_val[i,j]<-res1$estimate
              cor_p[i,j]<-res1$p.value
            }else{
              if(var.test(df.split[[1]]$c1, df.split[[2]]$c1)$p.value < 0.05){
                res1 <- t.test(df.split[[1]]$c1, df.split[[2]]$c1, var.equal = F)
                cor_val[i,j]<-abs(res1$estimate[1] - res1$estimate[2])
                cor_p[i,j]<-res1$p.value
              }
              else{
                res1 <- t.test(df.split[[1]]$c1, df.split[[2]]$c1, var.equal = T)
                cor_val[i,j]<-abs(res1$estimate[1] - res1$estimate[2])
                cor_p[i,j]<-res1$p.value
              }
            }
          }else{
            res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="spearman",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }

        }
        if(flag1 & flag2){
          res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="spearman",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }
      }
    }
  }
  
  if(method=="lm"){
    for(i in 1:ncol(df1)){
      for (j in 1:ncol(df2)) {
        
        var1 <- paste0("df2$",names(df2)[j])
        lm_variable <- names(df2)[-j]
        lm_variable <- paste0("df2$",lm_variable)
        lm_model <- paste(lm_variable,collapse = "+")
        lm_model <- paste0(var1 , "~","df1[,names(df1)[i]]" , "+",lm_model)
        try(res1<-lm(as.formula(lm_model) , na.action = na.omit), silent = TRUE)
        if(class(res1) != "try-error"){
          cor_val[i,j]<-as.numeric(lm.beta(res1)[1])
          cor_p[i,j]<-coef(summary(res1))[2,4]
        }
      }
    }
  }
  suppressMessages(library(reshape2))
  suppressWarnings(library(ggplot2))
  
  data <- cbind.data.frame(melt(cor_val),melt(cor_p)[,3])
  names(data) <- c("var1", "var2","cor_val","cor_p")
  if(Plot){
    #x_lab =names(df2)
    textMatrix = paste(signif(cor_val, 3), "\n(",
                       signif(cor_p, 3), ")", sep = "")
    dim(textMatrix) = dim(cor_val)
    
    data <- cbind.data.frame(data,melt(textMatrix)[,3])
    names(data)[5] = "text"
    
    cor_plot <- ggplot(data, aes(x = var2, y = var1, fill = cor_val)) +
      geom_tile(color = "black") +
      geom_text(aes(label = text), color = "black", size = 3) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",limits=c(-1,+1))+
      labs(x = "",y = "") +ggtitle(plot.title)+
      theme(legend.title = element_blank(), text = element_text(size = 20, color = "black"), plot.title = element_text(size = 14),
            legend.text = element_text(size = 12),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
  }
  if(return_melt){
    return(list(Result=data[,1:4],Plot = cor_plot))
  }else{
    colnames(cor_val) = paste0(colnames(cor_val) , ".corr")
    colnames(cor_p) = paste0(colnames(cor_p) , ".Pvalue")
    data <- cbind.data.frame(cor_val, cor_p)
    return(list(Result=data,Plot=cor_plot))
  }
  
}


scatter.smooth.with.lm <- function(x, y, xlabel="",ylabel="", title=""){
  library(ggplot2)
  data=cbind.data.frame(x,y)
  ggplot(data = data, aes_string(x = x, y = y)) + geom_point() + 
    geom_smooth(method = "lm", se = T) + ggtitle(title)+xlab(xlabel)+ylab(ylabel)
}
######################################################################################

setwd("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Brain EWAS Paper")

suppressMessages(library(WGCNA))
suppressMessages(library(stringr))
suppressMessages(library("QuantPsyc"))
suppressMessages(library(lsr))
suppressMessages(library(gridExtra))

net.file="WGCNA/March2024/Results/Pitts.Psycho.0.2.Pow3.wgcna.network.rds"
betas.file="BDR.data/BDR.All.Psycho.0.2.AD.braak.lmCorrected.betas.rds"
pheno.file="BDR.data/BDR.All.Psycho.0.2.AD.wgcna.pheno.csv"
trait="Psychosis"
#covars_fact <- "Sex,BraakStage,Plate,SentrixID,TissueType" 
covars_fact <- "Sex,BraakStage,Plate,SentrixID" 
covars_num <- "Age,CellProportion"
out_pref = "WGCNA/March2024/Results/BDR.0.2.lmCorrected.WGCNA.Module-Trait"
modules <- c("mediumpurple3","magenta1","firebrick3","lightcoral","lightcyan","darkgoldenrod1","green4","linen","mistyrose","seashell4","antiquewhite","moccasin","greenyellow","darkred","lightblue1")


analysis.type="tsppr"
# pr -> Pearson correlation
# sp -> Spearman correlation
# lm -> Linear regression
# prsp -> Pearson correlation for numeric variables and Spearman correlation for categorical variables
# tsppr -> Pearson correlation for two numeric variables,
#          Spearman correlation for two categorical variables or one categorical (more than one groups) and one numeric variables,
#          t-test for a binary and a numeric variable

calc_ME = T   # Calculate Module Eigengene from expression matrix
corr.plot= F  # correlation plot between MEs and Phenotype
scatter.plot = F # Scatter plot for MEs and Trait variable
save_csv=T # Save the results in csv format
all_moduls=T # Run the analysis on all modules
PCA_betas=F # PCA analysis and plot for expression matrix
###############################################################
pheno = read.csv(pheno.file,row.names = 1,stringsAsFactors = F)
betas = readRDS(betas.file)
net = readRDS(net.file)

betas = betas[rownames(betas) %in% names(net$colors),]
index = match(rownames(betas),names(net$colors))
net$colors = net$colors[index]
if(!identical(names(net$colors),rownames(betas))){
  print("CpG IDs in network and beta files are not match. It will be matched using shared CpGs...")
  index = intersect(names(net$colors),rownames(betas))
  betas = betas[index , ]
  net$colors = net$colors[index]
}
if(!identical(colnames(betas) , rownames(pheno))){
  print("Sample IDs in phenotype and beta files are not match. It will be matched using shared samples...")
  index = intersect(colnames(betas) , rownames(pheno))
  betas = betas[,index]
  pheno = pheno[index ,]
}

covars_fact=str_split(covars_fact,pattern = ',',simplify = T)[1,]
covars_num=str_split(covars_num,pattern = ',',simplify = T)[1,]

pheno.1 <- as.data.frame(as.numeric(as.factor(pheno[,trait])))

if(covars_fact[1] != "")
  for (i in 1:length(covars_fact)) {
    pheno.1 <- cbind.data.frame(pheno.1,as.numeric(as.factor(pheno[,covars_fact[i]])))
  }

if(covars_num[1] != "")
  for (i in 1:length(covars_num)) {
    pheno.1 <- cbind.data.frame(pheno.1,as.numeric(pheno[,covars_num[i]]))
  }

col.names = c(trait,covars_fact,covars_num)
col.names[col.names %in% ""] <- NA
names(pheno.1) <- col.names[!is.na(col.names)]

if(all_moduls){
  modules = unique(net$colors)
}

if(calc_ME){
  Mes = moduleEigengenes(expr = t(betas) , colors = net$colors,softPower = 3)$eigengenes
  net$MEs = Mes
}


cor_ <- matrix(data = NA, ncol = ncol(pheno.1), nrow = length(modules))
colnames(cor_) <- names(pheno.1)
rownames(cor_) <- modules
cor_pval <- cor_

colnames(net$MEs) <- str_remove(colnames(net$MEs),pattern = "ME")

result <- data.frame.relationship(df1 = net$MEs[,modules],df2 = pheno.1,method = analysis.type, 
                                  plot.title = paste("ME-Trait using",analysis.type,"correlation"),
                                  return_melt = F,Plot = corr.plot,categorical_columns = c(trait,covars_fact))
if(save_csv){
  write.csv(result$Result , file=paste0(out_pref,".",analysis.type,".csv"))
}
if(corr.plot){
  tiff(filename = paste0(out_pref,".",analysis.type,".tif"),units = "in", width = 10, height = 10,res = 500)
  print(result$Plot)
  graphics.off()
}
if(PCA_betas){
  pca = prcomp(t(betas), scale. = T , center = T)
  pca.betas = pca$x
  result <- data.frame.relationship(df1 = as.data.frame(pca.betas[,1:10]),df2 = pheno.1,method = analysis.type, 
                                    plot.title = paste(analysis.type,"correlation plot after removing batch effects"),
                                    return_melt = F,Plot = corr.plot,categorical_columns = c(trait,covars_fact))
  tiff(filename =  paste0(out_pref,".betas.",analysis.type,".tif"),units = "in", width = 10, height = 10,res = 500)
  print(result$Plot)
  graphics.off()

}

if(scatter.plot){
  p=list()
  for (i in 1:length(modules)) {
    p[[modules[i]]] = scatter.smooth.with.lm(pheno.1[,trait],net$MEs[,modules[i]],xlabel = "Psychosis",ylabel = "ME",title = modules[i])
  }
  
  tiff(filename = paste0(out_pref,".lm.scatter.tif"),units = "in", width = 10, height = 10,res = 500)
  grid.arrange(grobs=p)
  graphics.off()
}
