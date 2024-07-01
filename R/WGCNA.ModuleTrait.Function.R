data.frame.relationship <- function(df1 , df2, method="prsp", Plot=T,plot.title="",return_melt=T, categorical_columns){
  # method:
  #   p -> Pearson correlation 
  #   s -> Spearman correlation 
  #   l -> Linear regression 
  #   ps -> Pearson correlation for numeric variables and Spearman correlation for categorical variables 
  #   tps -> Pearson correlation for two numeric variables, Spearman correlation for two categorical variables 
  #          or one categorical (more than two groups) and one numeric variables, t-test for a binary and a numeric variable
  
  method <- match.arg(method,c("p","s","ps","l","tps"),several.ok = F)
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
  
  if(method=="p"){
    for(i in names(df1)){
      for (j in names(df2)) {
        res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="pearson",exact = FALSE)
        cor_val[i,j]<-res1$estimate
        cor_p[i,j]<-res1$p.value
      }
    }
  }
  if(method=="s"){
    for(i in names(df1)){
      for (j in names(df2)) {
        res1<-cor.test(as.numeric(df1[,i]),as.numeric(df2[,j]), method="spearman",exact = FALSE)
        cor_val[i,j]<-res1$estimate
        cor_p[i,j]<-res1$p.value
      }
    }
  }
  if(method=="ps"){
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
  if(method=="tps"){
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
  
  if(method=="l"){
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

