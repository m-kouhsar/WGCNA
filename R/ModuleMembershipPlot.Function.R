Module.Membership.Plot <- function(net.colors, expr.mat,metadata,trait,categorical.trait=TRUE,cofounders,modules,selected.genes=NA,
                                   pval.threshold=NA,p.adjust.threshodl=NA,p.adjust.method="bonferroni",soft.power,plot.title=""){
  suppressMessages(library(WGCNA))
  suppressMessages(library(ggplot2))
  suppressMessages(library(limma))
  
  if(!identical(names(net.colors),rownames(expr.mat))){
    warning("Genes/Probes in network colors and expression matrix are not identical. The intersection will be used.")
    index <- intersect(names(net.colors),rownames(expr.mat))
    net.colors <- net.colors[index]
    expr.mat <- expr.mat[index, ]
  }
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
    design <- model.matrix(as.formula(paste0("~",trait,"+",paste(cofounders,collapse = "+"))),data = metadata)
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
    message(paste("Generating Plot for",modules[i],"module..."))
    p <- ggplot(data = plot.data,aes(x=MM,y=GS))+
      geom_point(color="black")+
      geom_smooth(method = "lm") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))+
      xlab("Module Membership") +
      ylab("Probe Significance") 
    
    if(is.na(plot.title)){
      p <- p+ggtitle(paste0(modules[i],"\n","(Correlation=",round(c$estimate,2),", P-value=",formatC(c$p.value,format = "e",2),")"))
    }else{
      p <- p+ggtitle(paste0(plot.title,"\n","(Correlation=",round(c$estimate,2),", P-value=",formatC(c$p.value,format = "e",2),")"))
    }
    
    pval_na <- is.na(pval.threshold)
    padj_na <- is.na(p.adjust.threshodl)
    
    if(!pval_na & padj_na){
      selected.genes <- intersect(rownames())
    }else if(pval_na & !padj_na){
      
    }else if(!pval_na & !padj_na){
      
    }
    if(!is.na(selected.genes[1])){
      p <- p + geom_point(data=plot.data[plot.data$ID %in% selected.genes,],
                          pch=21, fill="blue", size=4, colour="blue", stroke=1)
    }
    Corrs$Module[i] <- modules[i]
    Corrs$Size[i] <- length(net.colors[net.colors==modules[i]])
    Corrs$Cor[i] <- c$estimate
    Corrs$Pvalue[i] <- c$p.value
    Plots[[i]] <- p
    Data[[i]] <- plot.data
  }
  return(list(Plots = Plots , Cor.test = Corrs , Data = Data))
}

