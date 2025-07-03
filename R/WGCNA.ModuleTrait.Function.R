ModuleTrait <- function(MEs , Pheno, method="cor", Plot=T,plot.title="",return_melt=T, 
                                    Trait , Factor_Covars, Numeric_Covars){
  # method:
  #        "lm" -> Linear regression adjusted for all covariates
  #        "cor" -> Pearson correlation for numeric variables and Spearman correlation for categorical variables
  #        "test" -> Pearson correlation for numeric variables
  #                 t-test for binary variables
  #                 anova test for categorical and numeric variables
  
  suppressMessages(library(QuantPsyc))
  method <- match.arg(method,c("cor","lm","test"),several.ok = F)
  
  if( !(Trait %in% colnames(Pheno))){
    stop("Trait variable ('", Trait,"') was not found in Pheno.")
  }
  if( !all(Factor_Covars %in% colnames(Pheno))){
    stop("The following Factor variables were not found in Pheno:\n",setdiff(Factor_Covars , colnames(Pheno)))
  }
  if( !all(Numeric_Covars %in% colnames(Pheno))){
    stop("The following Numerical variables were not found in Pheno:\n",setdiff(Numeric_Covars , colnames(Pheno)))
  }
  if(nrow(MEs) != nrow(Pheno)){
    stop("MEs and Pheno dataframes must have the same number of rows")
  }
  all_var <- c(Trait , Factor_Covars , Numeric_Covars)
  cor_val <- matrix(data = NA, nrow = ncol(MEs), ncol = length(all_var))
  colnames(cor_val) <- all_var
  rownames(cor_val) <- names(MEs)
  cor_p <- cor_val
  cor_plot <- NA
  
  if(method=="cor"){
    
    colnames(cor_val) <- paste0(colnames(cor_val),".Corr")
    
    for(i in 1:ncol(MEs)){
      
      for (j in 1:length(all_var)) {
        
        if(all_var[j] %in% c(Trait , Factor_Covars)){
          
          message("Spearman correlation test between ",colnames(MEs)[i] , " and ",all_var[j])
          res1<-cor.test(as.numeric(MEs[,i]),as.numeric(as.factor(Pheno[,all_var[j]])), method="spearman",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
          
        }else{
          message("Pearson correlation test between ",colnames(MEs)[i] , " and ",all_var[j])
          res1<-cor.test(as.numeric(MEs[,i]),as.numeric(Pheno[,all_var[j]]), method="pearson",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }
      }
    }
  }
  
  # https://stackoverflow.com/questions/52811684/running-a-two-sample-t-test-with-unequal-sample-size-in-r
  if(method=="test"){
    
    col.type <- vector(length = ncol(cor_val) , mode = "character")
    
    for(i in 1:ncol(MEs)){
      
      for (j in 1:length(all_var)) {
        
        if(all_var[j] %in% c(Trait , Factor_Covars)){
          
          df <- cbind.data.frame(module=MEs[,i],variable=Pheno[,all_var[j]])
          df.split <- split(df,as.factor(df$variable),drop = T)
          
          if(length(df.split) == 2){
            
            message("T test on ", colnames(MEs)[i], " ME grouped by ",all_var[j])
            col.type[j] = "t_stat"
            
            res1 <- t.test(df.split[[1]]$module, df.split[[2]]$module)
            cor_val[i,j]<-res1$statistic
            cor_p[i,j]<-res1$p.value
            
          }else{
            
            message("ANOVA test for ",colnames(MEs)[i] , " grouped by ",all_var[j])
            col.type[j] = "F_val"
            
            res1<-aov(formula = module~variable , data = df )
            cor_val[i,j]<-summary(res1)[[1]]["F value"][1,1]
            cor_p[i,j]<-summary(res1)[[1]]["Pr(>F)"][1,1]
            
          }
        }else{
          
          message("Pearson correlation test between ",colnames(MEs)[i] , " and ",all_var[j])
          col.type[j] = "Corr"
          
          res1<-cor.test(as.numeric(MEs[,i]),as.numeric(Pheno[,all_var[j]]), method="pearson",exact = FALSE)
          cor_val[i,j]<-res1$estimate
          cor_p[i,j]<-res1$p.value
        }
      }
    }
    
    colnames(cor_val) <- paste(colnames(cor_val) , col.type , sep = ".")
  }
  
  if(method=="lm"){
    
    cor_val <- vector(length = ncol(MEs) , mode = "list")
    cor_p <- vector(length = ncol(MEs) , mode = "list")
    names(cor_p) = names(cor_val) <- colnames(MEs)
    lm_data <- cbind.data.frame(MEs , Pheno)
    for(m in colnames(MEs)){
      
      message("Linear regression analaysis on ", m , " :")
      
      lm_model <- paste(m,paste(all_var , collapse = "+") , sep = "~")
      message("      lm model: ",lm_model)
      lm_model <- as.formula(lm_model)
      
      try(res1<-lm(lm_model ,data = lm_data, na.action = na.omit), silent = TRUE)
      if(class(res1) != "try-error"){
        lm_summary <- summary(res1)
        cor_val[[m]]<- lm_summary$coefficients[-1,1]
        cor_p[[m]]<- lm_summary$coefficients[-1,4]
      }else{
        message("lm failed!")
      }
    }
    cor_val <- as.matrix(do.call(rbind.data.frame, cor_val))
    cor_p <- as.matrix(do.call(rbind.data.frame, cor_p))
    colnames(cor_val) = colnames(cor_p) = names(lm_summary$coefficients[-1,1])
    rownames(cor_val) = rownames(cor_p) = colnames(MEs)
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
    colnames(cor_p) = paste0(colnames(cor_p) , ".Pval")
    data <- cbind.data.frame(cor_val, cor_p)
    return(list(Result=data,Plot=cor_plot))
  }
  
}
