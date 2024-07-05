Module.Trait <- function(ME , Pheno, outlier,method="prsp", Plot=T,plot.title="",return_melt=T, categorical_columns){
  # method:
  #   pr -> Pearson correlation
  #   sp -> Spearman correlation
  #   lm -> Linear regression
  #   prsp -> Pearson correlation for numeric variables and Spearman correlation for categorical variables
  #   anovapr -> Pearson correlation for numeric variables and one-way ANOVA test for caregorical variables
  
  method <- match.arg(method,c("pr","sp","prsp","lm","anovapr"),several.ok = F)
  if(!is.data.frame(ME)){
    stop("ME must be a dataframe")
  }
  if(!is.data.frame(ME)){
    stop("Pheno must be a dataframe")
  }
  
  cor_val <- matrix(data = NA, nrow = ncol(ME), ncol = ncol(Pheno))
  colnames(cor_val) <- names(Pheno)
  rownames(cor_val) <- names(ME)
  cor_p <- cor_val
  cor_plot <- NA
  
  if(method=="pr"){
    for(i in names(ME)){
      for (j in names(Pheno)) {
        if(!is.na(outlier[,i])){
          res1<-cor.test(as.numeric(ME[-outlier[,i],i]),as.numeric(Pheno[-outlier[,i],j]), method="pearson",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }else{
          res1<-cor.test(as.numeric(ME[,i]),as.numeric(Pheno[,j]), method="pearson",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }

      }
    }
  }
  if(method=="sp"){
    for(i in names(ME)){
      for (j in names(Pheno)) {
        if(!is.na(outlier[,i])){
          res1<-cor.test(as.numeric(ME[-outlier[,i],i]),as.numeric(Pheno[-outlier[,i],j]), method="spearman",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }else{
          res1<-cor.test(as.numeric(ME[,i]),as.numeric(Pheno[,j]), method="spearman",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }
      }
    }
  }
  if(method=="prsp"){
    for(i in names(ME)){
      for (j in names(Pheno)) {
        if((i %in% categorical_columns)|(j %in% categorical_columns)){
          if(!is.na(outlier[,i])){
            res1<-cor.test(as.numeric(ME[-outlier[,i],i]),as.numeric(Pheno[-outlier[,i],j]), method="spearman",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }else{
            res1<-cor.test(as.numeric(ME[,i]),as.numeric(Pheno[,j]), method="spearman",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }
        }else{
          if(!is.na(outlier[,i])){
            res1<-cor.test(as.numeric(ME[-outlier[,i],i]),as.numeric(Pheno[-outlier[,i],j]), method="pearson",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }else{
            res1<-cor.test(as.numeric(ME[,i]),as.numeric(Pheno[,j]), method="pearson",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }
        }
      }
    }
  }
  if(method=="anovapr"){
    for(i in names(ME)){
      for (j in names(Pheno)) {
        if((i %in% categorical_columns)|(j %in% categorical_columns)){
          if(!is.na(outlier[,i])){
            res1<-aov(ME[-outlier[,i],i] ~ Pheno[-outlier[,i],j])
            res2 <- as.data.frame(etaSquared(res1,anova = T))
            cor_val[i,j]<-res2$eta.sq[1]
            cor_p[i,j]<-res2$p[1]
          }else{
            res1<-aov(ME[,i] ~ Pheno[,j])
            res2 <- as.data.frame(etaSquared(res1,anova = T))
            cor_val[i,j]<-res2$eta.sq[1]
            cor_p[i,j]<-res2$p[1]
          }
        }else{
          if(!is.na(outlier[,i])){
            res1<-cor.test(as.numeric(ME[-outlier[,i],i]),as.numeric(Pheno[-outlier[,i],j]), method="pearson",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }else{
            res1<-cor.test(as.numeric(ME[,i]),as.numeric(Pheno[,j]), method="pearson",exact = FALSE)
            cor_val[i,j]<-res1$estimate
            cor_p[i,j]<-res1$p.value
          }
        }
      }
    }
  }
  if(method=="lm"){
    for(i in 1:ncol(ME)){
      for (j in 1:ncol(Pheno)) {
        if(!is.na(outlier[,i])){
          Pheno.1 = Pheno[-outlier[,i],]
          ME.1 = ME[-outlier[,i],]
          
          var1 <- paste0("Pheno.1$",names(Pheno.1)[j])
          lm_variable <- names(Pheno.1)[-j]
          lm_variable <- paste0("Pheno.1$",lm_variable)
          lm_model <- paste(lm_variable,collapse = "+")
          lm_model <- paste0(var1 , "~","ME.1[,names(ME.1)[i]]" , "+",lm_model)
          try(res1<-lm(as.formula(lm_model) , na.action = na.omit), silent = TRUE)
          if(class(res1) != "try-error"){
            cor_val[i,j]<-as.numeric(lm.beta(res1)[1])
            cor_p[i,j]<-coef(summary(res1))[2,4]
          }
        }else{
          var1 <- paste0("Pheno$",names(Pheno)[j])
          lm_variable <- names(Pheno)[-j]
          lm_variable <- paste0("Pheno$",lm_variable)
          lm_model <- paste(lm_variable,collapse = "+")
          lm_model <- paste0(var1 , "~","ME[,names(ME)[i]]" , "+",lm_model)
          try(res1<-lm(as.formula(lm_model) , na.action = na.omit), silent = TRUE)
          if(class(res1) != "try-error"){
            cor_val[i,j]<-as.numeric(lm.beta(res1)[1])
            cor_p[i,j]<-coef(summary(res1))[2,4]
          }
        }
      }
    }
  }
  suppressMessages(library(reshape2))
  suppressWarnings(library(ggplot2))
  
  data <- cbind.data.frame(melt(cor_val),melt(cor_p)[,3])
  names(data) <- c("var1", "var2","cor_val","cor_p")
  if(Plot){
    #x_lab =names(Pheno)
    textMatrix = paste(signif(cor_val, 3), "\n(",
                       signif(cor_p, 3), ")", sep = "")
    dim(textMatrix) = dim(cor_val)
    
    data <- cbind.data.frame(data,melt(textMatrix)[,3])
    names(data)[5] = "text"
    
    cor_plot <- ggplot(data, aes(x = var2, y = var1, fill = cor_val)) +
      geom_tile() +
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

net.file="WGCNA/mad.0.5.lm.v2.pow3/Pitts.All.Psycho.0.2.v2.lmCorrected.wgcna.network.rds"
betas.file="BDR.data/BDR.All.Psycho.0.2.AD.braak.lmCorrected.betas.rds"
pheno.file="BDR.data/BDR.All.Psycho.0.2.AD.pheno.csv"
outlier_file="WGCNA/mad.0.5.lm.v2.pow3/moduleOutlier.csv"
trait="Psychosis"
covars_fact <- "Sex,BraakStage,SentrixID,Plate" 
covars_num <- "Age,CellProportion"
out_pref = "Figures/wgcna.mad.0.5.v2.pow3/Module.Trait/BDR.0.2.lmCorrected"
#modules <- "mediumpurple3,magenta1,firebrick3,lightcoral,lightcyan,darkgoldenrod1,linen,lavenderblush,greenyellow,darkred,lightblue1"
modules <- "seashell4,mediumpurple3,thistle2,green4,firebrick3,sienna2,lightcoral,antiquewhite,magenta1,lightcyan,burlywood,linen,thistle3,lavenderblush,darkgoldenrod1,tan2,greenyellow,whitesmoke,orange,midnightblue,mistyrose,darkred,lightblue1"

analysis.type="sp"
# pr -> Pearson correlation
# sp -> Spearman correlation
# lm -> Linear regression
# prsp -> Pearson correlation for numeric variables and Spearman correlation for categorical variables
# anovapr -> Pearson correlation for numeric variables and one-way ANOVA test for caregorical variables
calc_ME = T
corr.plot= T
scatter.plot = T
save_csv=F
###############################################################
pheno = read.csv(pheno.file,row.names = 1,stringsAsFactors = F)
betas = readRDS(betas.file)
net = readRDS(net.file)
modules.outliers = read.csv(outlier_file, row.names = 1)

betas = betas[rownames(betas) %in% names(net$colors),]
index = match(rownames(betas),names(net$colors))
net$colors = net$colors[index]
if(!identical(names(net$colors),rownames(betas))){
  index = intersect(names(net$colors),rownames(betas))
  betas = betas[index , ]
  net$colors = net$colors[index]
}
if(!identical(colnames(betas) , rownames(pheno))){
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

modules = str_split(modules,pattern = ',',simplify = T)[1,]

if(calc_ME){
  Mes = moduleEigengenes(expr = t(betas) , colors = net$colors,softPower = 3)$eigengenes
  net$MEs = Mes
}


cor_ <- matrix(data = NA, ncol = ncol(pheno.1), nrow = length(modules))
colnames(cor_) <- names(pheno.1)
rownames(cor_) <- modules
cor_pval <- cor_

colnames(net$MEs) <- str_remove(colnames(net$MEs),pattern = "ME")

result <- Module.Trait(ME = net$MEs[,modules],Pheno = pheno.1,outlier = modules.outliers,method = analysis.type,
                       plot.title = paste("ME-Trait using",analysis.type,"correlation after removing module outliers"),
                                  return_melt = F,Plot = corr.plot,categorical_columns = c(trait,covars_fact))
if(save_csv){
  write.csv(result$Result , file=paste0(out_pref,".",analysis.type,".outlier.csv"))
}
if(corr.plot){
  tiff(filename = paste0(out_pref,".",analysis.type,".outlier.tif"),units = "in", width = 10, height = 10,res = 500)
  print(result$Plot)
  graphics.off()
}

if(scatter.plot){
  p=list()
  for (i in 1:length(modules)) {
    if(!is.na(modules.outliers[,i])){
      p[[modules[i]]] = scatter.smooth.with.lm(pheno.1[-modules.outliers[,i],trait],net$MEs[-modules.outliers[,i],modules[i]],xlabel = "Psychosis",ylabel = "ME",title = modules[i])
    }else{
      p[[modules[i]]] = scatter.smooth.with.lm(pheno.1[,trait],net$MEs[,modules[i]],xlabel = "Psychosis",ylabel = "ME",title = modules[i]) 
    }
  }
  
  tiff(filename = paste0(out_pref,".outlier.lm.scatter.tif"),units = "in", width = 10, height = 10,res = 500)
  grid.arrange(grobs=p)
  graphics.off()
}
