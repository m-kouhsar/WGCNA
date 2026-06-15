Module.Membership.Plot <- function(net.colors, expr.mat,metadata,trait,categorical.trait=TRUE,cofounders = NA,modules,
                                   legend_GS_pval=NA, label_GS_pval = NA , p.adjust.method="bonferroni",soft.power,plot.title=""){
  suppressMessages(library(WGCNA))
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggrepel))
  suppressMessages(library(limma))
  
  ##########################################################################################
  if(is.character(net.colors)){
    if(is.null(names(net.colors))){
      stop("Input Error: net.colors must be a named charachter vector.")
    }
  }else{
    stop("Input Error: net.colors must be a named charachter vector.")
  }
  
  if(is.na(legend_GS_pval)){
    legend_GS_pval <- 0
  }else{
    if (!is.numeric(legend_GS_pval) || legend_GS_pval < 0 || legend_GS_pval > 1) {
      stop("Input Error: 'legend_GS_pval' must be a numeric value between 0 and 1, or NA.")
    }
  }
  if(is.na(label_GS_pval)){
    label_GS_pval <- 0
  }else{
    if (!is.numeric(label_GS_pval) || label_GS_pval < 0 || label_GS_pval > 1) {
      stop("Input Error: 'label_GS_pval' must be a numeric value between 0 and 1, or NA.")
    }
  }
  
  if(!identical(names(net.colors),rownames(expr.mat))){
    warning("Genes/Probes in network colors and expression matrix are not identical. The intersection will be used.")
    index <- intersect(names(net.colors),rownames(expr.mat))
    if(length(index) <= 0){
      stop("There is no shared gene/probe between network colors and the expression matrix!")
    }else{
      net.colors <- net.colors[index]
      expr.mat <- expr.mat[index, ] 
    }
  }
  
  ####################################################################################################
  message("Calculating Modules Eigengenes")
  MEs = moduleEigengenes(expr = t(expr.mat) , colors = net.colors,softPower = soft.power)$eigengenes
  MEs <- orderMEs(MEs)
  names(MEs) = substring(names(MEs), 3)
  MEs <- MEs[,modules, drop = F]
  nSamples = ncol(expr.mat)
  
  message("Calculating Modules Memberships and P values...")
  geneModuleMembership = as.data.frame(WGCNA::cor(t(expr.mat), MEs, method = "pearson"))
  names(geneModuleMembership) <- modules
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  MMPvalue.adj = as.data.frame(matrix(NA, nrow = nrow(MMPvalue), ncol = ncol(MMPvalue)))
  rownames(MMPvalue.adj) <- rownames(MMPvalue)
  names(MMPvalue.adj) <- modules
  
  for (m in modules) {
    module_genes <- names(net.colors[net.colors == m])
    m_pval <- MMPvalue[module_genes , m]
    m_pval.adj <- p.adjust(m_pval , method = p.adjust.method)
    MMPvalue.adj[module_genes , m] <- m_pval.adj
  }
  
  
  message("Calculating Gene Significant and P values...")
  if(categorical.trait){
    
    metadata[,trait] <- as.factor(metadata[,trait])
    if(all(is.na(cofounders))){
      design <- model.matrix(as.formula(paste0("~",trait)),data = metadata)
    }else{
      design <- model.matrix(as.formula(paste0("~",trait,"+",paste(cofounders,collapse = "+"))),data = metadata)
    }
    
    fit <- lmFit(expr.mat, design)
    fit <- eBayes(fit)
    trait_column_name <- colnames(design)[2] 
    
    adjusted_betas <- as.data.frame(fit$coefficients)[, trait_column_name, drop = F]
    GS_raw <- abs(adjusted_betas)
    geneTraitSignificance<- (GS_raw - min(GS_raw, na.rm = T))/(max(GS_raw, na.rm = T) - min(GS_raw, na.rm = T))
    GSPvalue <- fit$p.value[, trait_column_name , drop = F]
    GSPvalue.adj = apply(GSPvalue,2,function(x){p.adjust(x,method = p.adjust.method)})
  }else{
    geneTraitSignificance = as.data.frame(WGCNA::cor(t(expr.mat), as.numeric(as.factor(metadata[,trait])), method = "pearson"))
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
    GSPvalue.adj = apply(GSPvalue,2,function(x){p.adjust(x,method = p.adjust.method)})
  }
  
  
  Plots <- vector(mode = "list",length = length(modules))
  names(Plots) <- modules
  Data <- vector(mode = "list",length = length(modules))
  names(Data) <- modules
  Corrs <- data.frame(Module=vector(mode = "character",length = length(modules)),
                      Size=vector(mode = "character",length = length(modules)),
                      Cor=vector(mode = "character",length = length(modules)),
                      Pvalue=vector(mode = "character",length = length(modules)))
  
  for (i in 1:length(modules)) {
    moduleGenes = (net.colors==modules[i])
    
    plot.data <- cbind.data.frame(ID = rownames(geneModuleMembership)[moduleGenes], 
                                  MM=abs(geneModuleMembership[moduleGenes , modules[i]]),
                                  MM.Pval = MMPvalue[moduleGenes , modules[i]],
                                  MM.Padj = MMPvalue.adj[moduleGenes , modules[i]],
                                  GS=abs(geneTraitSignificance[moduleGenes , 1]),
                                  GS.Pval = GSPvalue[moduleGenes , 1],
                                  GS.Padj = GSPvalue.adj[moduleGenes , 1])
    c <- cor.test(plot.data$MM,plot.data$GS,method = "pearson")
    
    plot.data$GS_sig <- "Not Significant"
    plot.data$GS_sig[plot.data$GS.Pval < legend_GS_pval] <- "Significant"
    
    message(paste("Generating Plot for",modules[i],"module..."))
    
    p <- ggplot()+geom_point(data = plot.data,aes(x=MM,y=GS,color = GS_sig))+
      scale_color_manual(name = "Related to trait", values = c("Not Significant" = "black" ,"Significant" = "chocolate2" ))+
      geom_smooth(data = plot.data,aes(x=MM,y=GS),method = "lm",colour = "darkgray") + 
      theme_minimal() + 
      theme(plot.title = element_text(hjust = 0.5))+
      xlab("Module Membership") +
      ylab("Gene Significance") 
    
    if(legend_GS_pval == 0){
      p <- p + theme(legend.position="none")
    }
    if(is.na(plot.title)){
      p <- p+ggtitle(paste0("Module Membership vs Gene Significance plot for ",modules[i]," baesd on ",trait),
                     subtitle = paste0("(Correlation=",round(c$estimate,2),", P-value=",formatC(c$p.value,format = "e",2),")"))
    }else{
      p <- p+ggtitle(plot.title,
                          subtitle = paste0("(Correlation=",round(c$estimate,2),", P-value=",formatC(c$p.value,format = "e",2),")"))
    }
    if(label_GS_pval == 0){
      selected.genes <- NA
    }else{
      selected.genes <- plot.data$ID[plot.data$GS.Pval < label_GS_pval]
    }
    if(!is.na(selected.genes[1])){
      p <- p + geom_point(data=plot.data[plot.data$ID %in% selected.genes,],
                          aes(x=MM,y=GS),
                          shape = 1,           # 1 = open circle
                          size = 4,            # Larger than base size (2.5) to act as a ring
                          color = "lightblue3",      # Color of the outline
                          stroke = 1.2,        # Thickness of the circle's line
                          show.legend = FALSE)+
        geom_label_repel(data=plot.data[plot.data$ID %in% selected.genes,],
                         aes(x=MM,y=GS , label = ID),
                         fill = "white",
                         color = "black",
                         box.padding = 0.5,
                         point.padding = 0.6,
                         show.legend = F,
                         max.overlaps = Inf,
                         force = 10)
    }
    Corrs$Module[i] <- modules[i]
    Corrs$Size[i] <- length(net.colors[net.colors==modules[i]])
    Corrs$Cor[i] <- c$estimate
    Corrs$Pvalue[i] <- c$p.value
    Plots[[i]] <- p
    Data[[i]] <- plot.data
  }
  Corrs <- Corrs[order(as.numeric(Corrs$Pvalue) , decreasing = F),]
  return(list(Plots = Plots , Cor.test = Corrs , Data = Data))
}

