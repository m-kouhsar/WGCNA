mahalanobis.outlier <- function(Data , method = "pca", plot.title=NA , seed = NA){
  
  suppressMessages(library(car))
  suppressMessages(library(dplyr))
  suppressMessages(library(Rtsne))
  if(!is.na(seed)){
    set.seed(seed = seed) 
  }
  method = match.arg(arg = method , choices = c("pca" , "tsne") , several.ok = F)
  
  if(method == "tsne"){
    tsne_out <- Rtsne(t(Data),dims = 2,)
    Data.2D <- data.frame(D1 = tsne_out$Y[,1], 
                          D2 = tsne_out$Y[,2])
    rownames(Data.2D) <- colnames(Data)
  }
  if(method == "pca"){
    pc <- prcomp(t(Data),scale. = T,center = T)
    pc.importance <- summary(pc)$importance[2,]
    Data.2D <- as.data.frame(pc$x[,1:2])
  }
  
  center_ <- colMeans(Data.2D)
  cov_ <- cov(Data.2D)
  
  Data.2D$mdist <- mahalanobis(
    x = Data.2D,
    center = center_,
    cov = cov_
  )
  
  cutoff <- qchisq(p = 0.95, df = 2)
  R <- sqrt(cutoff)
  
  ellipse_ <- car::ellipse(
    center = center_[1:2],
    shape = cov_[1:2,1:2],
    radius = R,
    segments = 150,
    draw = FALSE
  )
  ellipse_ <- as.data.frame(ellipse_)
  colnames(ellipse_) <- colnames(Data.2D)[1:2]
  
  Data.2D$p <- pchisq(Data.2D$mdist, df = 2, lower.tail = FALSE)
  
  Data.2D <- Data.2D %>%
    mutate(Outlier = ifelse(mdist > cutoff, 'Yes', 'No'))
  if(method == "pca"){
    p <- ggplot(Data.2D, aes(x = PC1 , y = PC2, color = Outlier))
    plot.subtitle = "Dimentinality reduction method: PCA"
  }
  if(method == "tsne"){
    p <- ggplot(Data.2D, aes(x = D1 , y = D2, color = Outlier))
    plot.subtitle = "Dimentinality reduction method: tSNE"
  }
  if(is.na(plot.title)){
    plot.title = "Outlier detection using Mahalanobis distance and Chi-Squared Distribution"
  }
  p <- p +
    geom_point(size = 3) +
    geom_point(aes(center_[1], center_[2]) , size = 5 , color = 'blue') +
    geom_polygon(data = ellipse_, fill = 'white', color = 'black', alpha = 0.3) +
    scale_color_manual(values = c('gray44', 'red')) +
    labs(title = plot.title,subtitle = plot.subtitle) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list(Data.2D = Data.2D , Plot=p))
}
